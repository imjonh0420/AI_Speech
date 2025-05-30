/*
 * AI_Speech_Spectral_Envelope_Dynamic_Lifter.c
 * Author : Imjonh (Modified by AI Assistant)
 * Date : 2025-05-30
 * Description : Find the Spectral Envelope of Voice Signals based on the Cepstrum
 * using dynamic lifter points (1 to 319) and outputting all
 * envelopes to a single file.
 * Homework #3 for "AI and Speech Signal Processing"
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h> // For fabs

#define SAMPLE_SIZE 320
#define PI 3.14159265358979323846
#define EPSILON 1e-12 // To avoid log(0)

// DFT 함수
// arguments: DFT 할 신호, 실수부 출력 배열, 허수부 출력 배열, 샘플 크기 N
void Dft(float* signal, float* real_out, float* imag_out, int N) {
    for (int k = 0; k < N; k++) {
        real_out[k] = 0;
        imag_out[k] = 0;
        for (int n = 0; n < N; n++) {
            double angle = 2 * PI * k * n / N;
            real_out[k] += signal[n] * cos(angle);
            imag_out[k] -= signal[n] * sin(angle);
        }
    }
}

// IDFT 함수
void Idft(float* real_in, float* imag_in, float* signal_out, int N) {
    for (int n = 0; n < N; n++) {
        signal_out[n] = 0;
        for (int k = 0; k < N; k++) {
            double angle = 2 * PI * k * n / N;
            signal_out[n] += real_in[k] * cos(angle) - imag_in[k] * sin(angle);
        }
        signal_out[n] /= N;
    }
}

int main() {

    FILE* male_voice_data = NULL;
    FILE* output_all_envelopes_file = NULL;

    // ==== [ 정적 배열 선언 및 초기화 ] ====
    // 입력 받을 음성 신호 저장용 변수
    float speech_signal[SAMPLE_SIZE] = { 0.0f, };
    short data_short = 0;

    // 스펙트럼 관련 변수들
    float spectrum_real[SAMPLE_SIZE] = { 0.0f, };
    float spectrum_imag[SAMPLE_SIZE] = { 0.0f, };
    float spectrum_mag[SAMPLE_SIZE] = { 0.0f, };
    float log_spectrum_mag[SAMPLE_SIZE] = { 0.0f, };

    // cepstrum 관련 변수들
    float cepstrum_imag_input[SAMPLE_SIZE] = { 0.0f, }; // IDFT의 허수부 입력은 0
    float cepstrum[SAMPLE_SIZE] = { 0.0f, };
    float liftered_cepstrum[SAMPLE_SIZE] = { 0.0f, };

    // envelope 관련 변수들
    float envelope_real[SAMPLE_SIZE] = { 0.0f, };
    float envelope_imag[SAMPLE_SIZE] = { 0.0f, };
    float log_spectral_envelope[SAMPLE_SIZE] = { 0.0f, };

    // 파일 열기
    if (fopen_s(&male_voice_data, "Male.raw", "rb") != 0 || male_voice_data == NULL) {
        printf("Error opening Male.raw\n");
        return 1;
    }

    if (fopen_s(&output_all_envelopes_file, "spectral_envelope_all_lifter_points.txt", "w") != 0 || output_all_envelopes_file == NULL) {
        printf("Error opening spectral_envelope_all_lifter_points.txt\n");
        if (male_voice_data) fclose(male_voice_data);
        return 1;
    }

    // 1. 음성 신호 읽기 (첫 320 샘플)
    for (int n = 0; n < SAMPLE_SIZE; n++) {
        if (fread(&data_short, sizeof(short), 1, male_voice_data) != 1) {
            printf("Error reading from Male.raw or EOF reached prematurely at sample %d\n", n);
            fclose(male_voice_data);
            fclose(output_all_envelopes_file);
            return 1;
        }
        speech_signal[n] = (float)data_short;
    }
    fclose(male_voice_data);

    // 2. Hamming Window 적용
    for (int n = 0; n < SAMPLE_SIZE; n++) {
        speech_signal[n] *= (0.54f - 0.46f * cos(2.0f * PI * n / (SAMPLE_SIZE - 1.0f)));
    }

    // 3. 입력 신호 Male_voice에 대해 DFT를 수행하여 스펙트럼 S(k) 계산
    Dft(speech_signal, spectrum_real, spectrum_imag, SAMPLE_SIZE);

    // 4. Cepstrum 분석을 위해 전체 Sample Size에 대해 Spectrum magnitude 계산
    for (int k = 0; k < SAMPLE_SIZE; k++) {
        spectrum_mag[k] = sqrt( (spectrum_real[k] * spectrum_real[k]) + (spectrum_imag[k] * spectrum_imag[k]) );
        log_spectrum_mag[k] = log(spectrum_mag[k] + EPSILON);
    }

    // 5. 덧셈 연산을 위해 log를 도입한 spectrum magnitude를 IDFT하여 Low-time Liftering 준비 (켑스트럼 계산)
    // log_spectrum_mag는 실수이므로, IDFT의 입력에서 허수부는 0으로 설정 (cepstrum_imag_input)
    Idft(log_spectrum_mag, cepstrum_imag_input, cepstrum, SAMPLE_SIZE);

    // 리프터 포인트를 1부터 319까지 변경하면서 엔벨로프 계산 및 출력
    for (int current_lifter_point = 1; current_lifter_point < SAMPLE_SIZE; current_lifter_point++) {
        // liftered_cepstrum 배열을 0.0f로 초기화
        for(int i = 0; i < SAMPLE_SIZE; i++) {
            liftered_cepstrum[i] = 0.0f;
        }

        // 6. Low-time liftering, 대칭성 유지
        liftered_cepstrum[0] = cepstrum[0]; // DC 성분(0번째 켑스트럼 계수)은 항상 포함
        for (int n = 1; n < current_lifter_point; n++) {
            // current_lifter_point가 1이면 이 루프는 실행되지 않음.
            // n은 1부터 current_lifter_point-1 까지만큼의 계수를 가져옴.
            if (n < SAMPLE_SIZE) { // 양의 quefrency 축 계수
                liftered_cepstrum[n] = cepstrum[n];
            }
            if ((SAMPLE_SIZE - n) > 0 && (SAMPLE_SIZE - n) < SAMPLE_SIZE) { // 음의 quefrency 축 계수 (대칭)
                                                                         // (SAMPLE_SIZE - n)이 0이 되지 않도록 (이미 liftered_cepstrum[0] 처리됨)
                liftered_cepstrum[SAMPLE_SIZE - n] = cepstrum[SAMPLE_SIZE - n];
            }
        }

        // 7. log|H(k)|를 얻기 위해 liftered_cepstrum에 대해 DFT 수행
        Dft(liftered_cepstrum, envelope_real, envelope_imag, SAMPLE_SIZE);

        // 8. 게인 맞추기 (DC 성분 일치)
        float dc_offset_adjustment = 0;
        if (SAMPLE_SIZE > 0 ) { // envelope_real[0]이 0에 가까울 경우도 고려
            dc_offset_adjustment = log_spectrum_mag[0] - envelope_real[0];
        }

        // 8-1. log_spectral_envelope 게인 맞추고 파일 출력 (요청하신 160개만)
        for (int k = 0; k < SAMPLE_SIZE / 2; k++) {
            log_spectral_envelope[k] = envelope_real[k] + dc_offset_adjustment;
            fprintf(output_all_envelopes_file, "%f\n", log_spectral_envelope[k]);
        }
    } // End of current_lifter_point loop

    fclose(output_all_envelopes_file);

    printf("Processing complete. Output file generated:\n");
    printf(" - spectral_envelope_all_lifter_points.txt (All log-magnitude spectral envelopes for lifter points 1 to 319)\n");

    return 0;
}
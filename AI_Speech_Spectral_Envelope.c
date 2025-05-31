#include <stdio.h>
#include <math.h>

#define SAMPLE_SIZE 320
#define PI 3.14159265358979323846
#define EPSILON 1e-12 // To avoid log(0)

// DFT 함수
// arguments: DFT 할 신호, 실수부 출력 배열 주소, 허수부 출력 배열 주소, 샘플 크기 N
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
// arguments: 실수부 입력 배열, 허수부 입력 배열, IDFT 결과 출력 배열, 샘플 크기 N
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

    // 파일 포인터 선언, 초기화화
    FILE* male_voice_data = NULL;
    FILE* output_all_envelopes_file = NULL;

    // ==== [ 배열 선언 및 초기화 ] ====
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

    if (fopen_s(&output_all_envelopes_file, "spectral_envelope_all_liftering_points.txt", "w") != 0 || output_all_envelopes_file == NULL) {
        printf("Error opening spectral_envelope_all_liftering_points.txt\n");
        if (male_voice_data) fclose(male_voice_data);
        return 1;
    }

    // [1]. 음성 신호 읽기 (첫 320 샘플)
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

    // [1-1]. Hamming Window 적용
    for (int n = 0; n < SAMPLE_SIZE; n++) {
        double angle = 2.0 * PI * n / (SAMPLE_SIZE - 1.0f);
        speech_signal[n] *= (0.54f - 0.46f * cos(2.0f * PI * n / (SAMPLE_SIZE - 1.0f)));  // Hamming Window
        //speech_signal[n] *= 0.5f - 0.5f * cos(angle); // Hanning Window
        //speech_signal[n] *= (0.42f - 0.5f * cos(angle) + 0.08f * cos(2.0f * angle)); // Blackman Window
        //speech_signal[n] *= sin((PI/(SAMPLE_SIZE - 1.0))*(n+0.5)); // Sine Window
    }

    // [2]. 입력 신호 Male_voice에 대해 DFT를 수행하여 스펙트럼 S(k) 계산
    Dft(speech_signal, spectrum_real, spectrum_imag, SAMPLE_SIZE);

    // [3]. Cepstrum 분석을 위해 전체 Sample Size에 대해 Spectrum magnitude 계산
    for (int k = 0; k < SAMPLE_SIZE; k++) {
        spectrum_mag[k] = sqrtf( (spectrum_real[k] * spectrum_real[k]) + (spectrum_imag[k] * spectrum_imag[k]) ); // sqrtf 사용
        log_spectrum_mag[k] = logf(spectrum_mag[k] + EPSILON); // logf 사용
    }

    // [4]. 덧셈 연산을 위해 log를 도입한 spectrum magnitude를 IDFT하여 Low-time Liftering 준비 (켑스트럼 계산)
    Idft(log_spectrum_mag, cepstrum_imag_input, cepstrum, SAMPLE_SIZE);

    // [5]. 리프터링 포인트를 1부터 SAMPLE_SIZE-1까지 변경하면서 Low-time Liftering 수행 (엔벨로프 계산)
    for (int current_liftering_point = 1; current_liftering_point < SAMPLE_SIZE; current_liftering_point++) {
        // 5-1. liftered_cepstrum 배열을 0.0f로 초기화 (매 반복마다 수행)
        // cepstrum에서 lifering Point 이후를 0으로 초기화 하는 것이 아니라,
        // liftered_cepstrum을 0으로 초기화 해놓고 liftering Point 이전의 계수만 복사하는 방식
        for(int i = 0; i < SAMPLE_SIZE; i++) {
            liftered_cepstrum[i] = 0.0f;
        }

        // [5-2]. Low-time liftering, 대칭성 유지
        liftered_cepstrum[0] = cepstrum[0]; // 0번째 켑스트럼 계수는는 항상 포함됨.
        // current_liftering_point가 1이면 아래 루프는 실행되지 않음 (n=1부터 시작하고 n < 1 조건이 false)
        // n은 1부터 current_liftering_point-1 까지의 계수를 복사
        // cepstrum[SAMPLE_SIZE-n]은 cepstrum[n]과 같음 (대칭성)
        for (int n = 1; n < current_liftering_point; n++) {
            liftered_cepstrum[n] = cepstrum[n];
            liftered_cepstrum[SAMPLE_SIZE - n] = cepstrum[SAMPLE_SIZE - n];
        }

        // [6]. log|H(k)|를 얻기 위해 liftered_cepstrum에 대해 DFT 수행
        Dft(liftered_cepstrum, envelope_real, envelope_imag, SAMPLE_SIZE);

        // [6-1]. 게인 맞추기 (첫 번째 샘플 값을 일치 시킴.)
        // 이후, Excel을 이용해 추가적으로 Offset을 맞추는 작업을 수행 함
        float dc_offset_adjustment = 0;
        dc_offset_adjustment = log_spectrum_mag[0] - envelope_real[0];
        
        // [6-2]. log_spectral_envelope 게인 맞추고 파일 출력 (SAMPLE_SIZE / 2 개 포인트)
        for (int k = 0; k < SAMPLE_SIZE / 2; k++) {
            log_spectral_envelope[k] = envelope_real[k] + dc_offset_adjustment;
            fprintf(output_all_envelopes_file, "%f\n", log_spectral_envelope[k]);
        }
    } // End of current_liftering_point loop

    fclose(output_all_envelopes_file);
    printf("Processing complete. Output file generated:\n");
    printf(" - spectral_envelope_all_liftering_points.txt (All log-magnitude spectral envelopes for lifter points 1 to %d)\n", SAMPLE_SIZE - 1);

    return 0;
}
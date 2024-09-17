%% Part 4: SSS search & Cell Id Detection

FFT_size = 4096;
CP_length = 288;
SCS = 30e3;
Ts = 1/FFT_size/SCS;
CP_OFDM_length = FFT_size+CP_length;
num_sc = 240;
N_id_1 = 77;
N_id_2 = 2;

%% OFDM Modulation
SSS_stream = SSS_BPSK(N_id_1,N_id_2);                                             

% Map symbol to subcarrier
d_SSS = [zeros(56,1);SSS_stream;zeros(FFT_size-183,1)];
% FFT
OFDM_SSS_body = ifft(d_SSS)*sqrt(FFT_size);
% Add CP
CP_OFDM_SSS = [OFDM_SSS_body(end-CP_length+1:end);OFDM_SSS_body];

%% Channel and Noise
h = 1;
signal_after_channel = conv(CP_OFDM_SSS,h);

SNR = 20;
N_0 = 10^(-SNR/10) * (norm(signal_after_channel)^2/length(signal_after_channel));
noise = sqrt(N_0/2)*(randn(length(signal_after_channel),1) + 1j*randn(length(signal_after_channel),1));
received_SSS_signal = signal_after_channel + noise;

%% SSS search
corr = zeros(1,336);

for i = 0:335
    SSS_ref_stream = SSS_BPSK(i,N_id_2);
    d_SSS_ref = [zeros(56,1);SSS_ref_stream;zeros(FFT_size-183,1)];
    OFDM_SSS_ref_body = ifft(d_SSS_ref);
    CP_OFDM_SSS_ref = [OFDM_SSS_ref_body(end-CP_length+1:end);OFDM_SSS_ref_body];
    corr(i+1) = abs(received_SSS_signal' * CP_OFDM_SSS_ref) ;
end
[~,N_id_1_est_pos] = max(corr);
N_id_1_est = N_id_1_est_pos - 1;

figure;
plot(0:335,corr);
xlabel("N^1_{id} candidates")
ylabel("Correlation")


function BPSK_stream = SSS_BPSK(N_id_1,N_id_2)
    x_0 = zeros(127,1);
    x_1 = zeros(127,1);
    BPSK_stream = zeros(127,1);
    x_init = [1 0 0 0 0 0 0];
    x_0(1:7) = x_init;
    x_1(1:7) = x_init;
    
    for i = 1:120
        x_0(i+7) = mod(x_0(i+4)+x_0(i),2);
        x_1(i+7) = mod(x_1(i+1)+x_1(i),2);
    end
    for n = 0:126
        m_0 = mod(n + 15* floor(N_id_1/112) + 5*N_id_2,127);
        m_1 = mod(n + mod(N_id_1,112),127);
        BPSK_stream(n+1) = (1-2*x_0(m_0+1))*(1-2*x_1(m_1+1));
    end
end
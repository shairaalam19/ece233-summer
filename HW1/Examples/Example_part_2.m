%% Part 2: Frequency offset compensation & PSS search

FFT_size = 4096;
CP_length = 288;
SCS = 30e3;
Ts = 1/FFT_size/SCS;
CP_OFDM_length = FFT_size+CP_length;
num_sc = 240;
N_id_2 = 2;
delta_f = 2.03e3; % kHz

%% OFDM Modulation
PSS_stream = PSS_BPSK(N_id_2);                                             

% Map symbol to subcarrier
d_PSS = [zeros(56,1);PSS_stream;zeros(FFT_size-183,1)];
% FFT
OFDM_PSS_body = ifft(d_PSS)*sqrt(FFT_size);
% Add CP
CP_OFDM_PSS = [OFDM_PSS_body(end-CP_length+1:end);OFDM_PSS_body];

%% Channel and Noise
h = 1;
signal_after_channel = conv(CP_OFDM_PSS,h);

SNR = 20;
N_0 = 10^(-SNR/10) * (norm(signal_after_channel)^2/length(signal_after_channel));
noise = sqrt(N_0/2)*(randn(length(signal_after_channel),1) + 1j*randn(length(signal_after_channel),1));
received_PSS_signal = signal_after_channel .* exp(1j*2*pi*delta_f*(0:length(signal_after_channel)-1)'*Ts) + noise;

%% PSS search + coarse frequency offset estimation
freq_offset = -SCS/2:100:SCS/2;
corr = zeros(3,length(freq_offset));

for i = 0:2
    PSS_ref_stream = PSS_BPSK(i);
    d_PSS_ref = [zeros(56,1);PSS_ref_stream;zeros(FFT_size-183,1)];
    OFDM_PSS_ref_body = ifft(d_PSS_ref);
    CP_OFDM_PSS_ref = [OFDM_PSS_ref_body(end-CP_length+1:end);OFDM_PSS_ref_body];
    for j = 1:numel(freq_offset)
        corr(i+1,j) = abs(received_PSS_signal' * (CP_OFDM_PSS_ref .* ...
        exp(1j*2*pi*freq_offset(j)*(0:length(CP_OFDM_PSS_ref)-1)'*Ts)));
    end
end
[~,N_id_2_est_pos] = max(max(abs(corr),[],2));
N_id_2_est = N_id_2_est_pos - 1;

figure;
plot(freq_offset/1e3,corr(1,:));
hold on;
plot(freq_offset/1e3,corr(2,:));
plot(freq_offset/1e3,corr(3,:));
legend("N^2_{id} = 0","N^2_{id} = 1","N^2_{id} = 2")
xlabel("coarse frequency offset candidates (kHz)")
ylabel("Correlation")

%% Fine Frequency offset estimation
freq_offset_est = -angle(received_PSS_signal(FFT_size+1:FFT_size+CP_length)' ...
                  * received_PSS_signal(1:CP_length))/2/pi*SCS;



function BPSK_stream = PSS_BPSK(N_id_2)
    x = zeros(127,1);
    BPSK_stream = zeros(127,1);
    x_init = [0 1 1 0 1 1 1];
    x(1:7) = x_init;
    for i = 1:120
        x(i+7) = mod(x(i+4)+x(i),2);
    end
    for n = 0:126
        m = mod(n + 43*N_id_2,127);
        BPSK_stream(n+1) = 1-2*x(m+1);
    end
end
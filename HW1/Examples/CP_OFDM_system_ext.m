%% CP_OFDM_system.m
FFT_size = 4096;
CP_length = 288;
CP_OFDM_length = FFT_size+CP_length;
num_sc = 240;
SNR = -10:2:10;
num_MC = 1000;

%% OFDM Modulation
% QPSK symbol
QAM_mod = 4;
pilot_bit_stream = randi([0 1],num_sc*log2(QAM_mod),1);
QPSK_pilot_stream = qammod(pilot_bit_stream,QAM_mod,InputType='bit',UnitAveragePower=true);

data_bit_stream = randi([0 1],num_sc*log2(QAM_mod),1);                                             
QPSK_stream = qammod(data_bit_stream,QAM_mod,InputType='bit',UnitAveragePower=true);

% Map QPSK symbol to subcarrier
pilot = [QPSK_pilot_stream;zeros(FFT_size-num_sc,1)];
data = [QPSK_stream;zeros(FFT_size-num_sc,1)];
% FFT
OFDM_pilot_body = ifft(pilot);
OFDM_data_body = ifft(data);
% Add CP
CP_OFDM_pilot = [OFDM_pilot_body(end-CP_length+1:end);OFDM_pilot_body];
CP_OFDM_data = [OFDM_data_body(end-CP_length+1:end);OFDM_data_body];
% Concatenation
CP_OFDM_chain = [CP_OFDM_pilot;CP_OFDM_data];

figure;
plot(abs(CP_OFDM_chain));

% Channel
h = [1 0.5]';
received_signal = conv(CP_OFDM_chain,h);

% Monte-Carlo simulation
num_error_bit = zeros(num_MC,length(SNR));
for SNR_id = 1:length(SNR)
    for MC_id = 1:num_MC 
        % Noise
        received_signal_noisy = awgn(received_signal,SNR(SNR_id),'measured');
        
        % OFDM Demodulation
        received_CP_OFDM_chain = received_signal_noisy(1:CP_OFDM_length*2);
        
        % Remove CP
        received_OFDM_pilot_body = received_CP_OFDM_chain(CP_length+1:CP_OFDM_length);
        received_OFDM_data_body = received_CP_OFDM_chain(CP_OFDM_length+CP_length+1:end);
        
        % FFT
        received_pilot = fft(received_OFDM_pilot_body);
        received_data = fft(received_OFDM_data_body);
        
        % Channel Estimation & Equalization
        % ZF channel estimation
        ch_est = received_pilot(1:num_sc)./QPSK_pilot_stream;
        % Equalization
        equalized_data = received_data(1:num_sc)./ch_est;
        
        % QPSK demodulation
        demod_bit_stream = qamdemod(equalized_data,QAM_mod,OutputType='bit',UnitAveragePower=true);
        
        % Bit error
        num_error_bit(MC_id,SNR_id) = sum(demod_bit_stream~=data_bit_stream);   
    end
end
BER = sum(num_error_bit,1)/num_MC/num_sc/log2(QAM_mod);

% Plot
figure;
semilogy(SNR,BER)
grid on;
xlabel("SNR (dB)")
ylabel("BER")


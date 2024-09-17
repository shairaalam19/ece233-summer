%% CP_OFDM_system.m
FFT_size = 4096;
CP_length = 288;
CP_OFDM_length = FFT_size+CP_length;
num_sc = 240;

%% OFDM Modulation
% QPSK symbol
QAM_mod = 4;
pilot_bit_stream = randi([0 1],num_sc*log2(QAM_mod),1);
QPSK_pilot_stream = qammod(pilot_bit_stream,QAM_mod,InputType='bit',UnitAveragePower=true);

data_bit_stream = randi([0 1],num_sc*log2(QAM_mod),1);                                             
QPSK_stream = qammod(data_bit_stream,QAM_mod,InputType='bit',UnitAveragePower=true);

% Map symbol to subcarrier
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

%% Channel
h = [1 0.5]';
received_signal = conv(CP_OFDM_chain,h);

%% OFDM Demodulation
received_CP_OFDM_chain = received_signal(1:CP_OFDM_length*2);

% Remove CP
received_OFDM_pilot_body = received_CP_OFDM_chain(CP_length+1:CP_OFDM_length);
received_OFDM_data_body = received_CP_OFDM_chain(CP_OFDM_length+CP_length+1:end);

% FFT
received_pilot = fft(received_OFDM_pilot_body);
received_data = fft(received_OFDM_data_body);

%% Channel Estimation & Equalization
% ZF channel estimation
ch_est = received_pilot(1:num_sc)./QPSK_pilot_stream;
% Equalization
equalized_data = received_data(1:num_sc)./ch_est;

% QPSK demodulation
demod_bit_stream = qamdemod(equalized_data,QAM_mod,OutputType='bit',UnitAveragePower=true);


isequal(demod_bit_stream,data_bit_stream)

figure;
plot(real(equalized_data(1:num_sc)),imag(equalized_data(1:num_sc)),'.')
xlabel("I-channel")
ylabel("Q-channel")



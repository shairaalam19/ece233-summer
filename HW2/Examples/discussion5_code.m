%% Discussion 5) P1: Beamforming - Part a
clear all
N_T = [4, 16, 32];        %Number of transmit antennas
f_c = 2.4e9;              %Carrier frequency (Hz)
lambda = 3e8/f_c;         %Wavelength
d = lambda/2;             %Antenna spacing
phi = 0;                  %Steering angle (radians)
theta = -2*pi:0.001:2*pi;  %Angles that we will use to compute beamforming gains

figure
for n=1:length(N_T)
    gain = zeros(1,length(theta));                                                      %Initialize gain vector   
    w_BF = exp(-1j*2*pi*(0:N_T(n)-1)'*d*sin(phi)/lambda);                               %Beamforming vector
    for t = 1:length(theta)
        a_theta = exp(-1j*2*pi*(0:N_T(n)-1)'*d*sin(theta(t))/lambda)/sqrt(N_T(n));      %ULA array response  
        gain(t) = abs(w_BF'*a_theta)^2;                                                 %Array gain
    end
    polarplot(theta, gain, 'LineWidth', 1)
    hold on
end
legend('N_T=4', 'N_T=16', 'N_T=32')
rlim('auto')
grid on
thetalim([-90 90])
%--------------------------------------------------------------------------
    
%% Discussion 5) P1: Beamforming - Part b
clear all
N_T = 1:32;                     %Number of transmit antennas
f_c = 2.4e9;                    %Carrier frequency (Hz)
lambda = 3e8/f_c;               %Wavelength
d = lambda/2;                   %Antenna spacing
phi = pi/4;                     %Steering angle (radians)
theta = phi;                    %Angles that we will use to compute beamforming gains
gain = zeros(1,length(N_T));    %Initialize gain vector
for n = 1:length(N_T)
    w_BF = exp(-1j*2*pi*(0:N_T(n)-1)'*d*sin(phi)/lambda);                        %Beamforming vector
    a_theta = exp(-1j*2*pi*(0:N_T(n)-1)'*d*sin(theta)/lambda)/sqrt(N_T(n));      %ULA array response
    gain(n) = abs(w_BF'*a_theta)^2;                                              %Array gain
end

figure
plot(N_T, 10*log10(gain), N_T, 10*log10(N_T), 'o', 'LineWidth', 1)
xlabel('Number of antennas')
ylabel('Array gain (dB)')
grid on
legend('Array gain', '10log_{10}(N_T)')
%--------------------------------------------------------------------------

%% Discussion 5) P1: Beamforming - Part c
clear all
N_T = [4, 32];                                  %Number of transmit antennas
f_c = 2.4e9;                                    %Carrier frequency (Hz)
lambda = 3e8/f_c;                               %Wavelength
d = lambda/2;                                   %Antenna spacing
phi = pi/4;                                     %Steering angle (radians)
theta = phi;                                    %Angles that we will use to compute beamforming gains
sigma_e = 0:0.01:pi/8;                          %Standard deviation for the error in steering vector
iter = 10000;                                   %Number of iterations
gain = zeros(length(N_T), length(sigma_e));     %Initialize gain vector
for n = 1:length(N_T)
    a_theta = exp(-1j*2*pi*(0:N_T(n)-1)'*d*sin(theta)/lambda)/sqrt(N_T(n));   %ULA array response
    for k = 1:length(sigma_e)
        iterGain = zeros(1,iter);
        for t = 1:iter
            w_BF = exp(-1j*2*pi*(0:N_T(n)-1)'*d*sin(phi+sigma_e(k)*randn)/lambda);      %Beamforming vector
            iterGain(t)= abs(w_BF'*a_theta)^2;                                          %Array gain
        end
        gain(n,k) = mean(iterGain);
    end
end

figure
plot(sigma_e*180/pi, 10*log10(gain))
xlabel('\sigma_\epsilon (degrees)')
ylabel('Array gain (dB)')
legend('N_T=4', 'N_T=32')
%--------------------------------------------------------------------------



%% Discussion 5) P2: Multiplexing
clear all
N_T = 2;                                                            %No. of transmit antennas
N_R = 2;                                                            %No. of receive antennas 
M = 4;                                                              %Modulation level
Eb_N0 = 0:20;                                                       %Range of EB/N_0 
[avgBER_ZF, avgBER_MMSE, avgBER_SVD] = deal(zeros(size(Eb_N0)));    %Initialize
iter = 1000;                                                        %Number of iterations
%Loop over every value of Eb/No
for i = 1:length(Eb_N0)
    i
    %Convert Eb/No to linear scale
    Eb_N0_lin = 10^(Eb_N0(i)/10);
    %Convert Eb/No to SNR
    SNR_lin = 2*Eb_N0_lin; % for 4-QAM Es = 2Eb and Es/N0 = SNR if we assume raised cosine pulse
    for t = 1:iter
        %Generate random symbols
        b_in = randi([0 1],N_T*log2(M),1);
        x = qammod(b_in, M, 'InputType', 'bit', 'UnitAveragePower', true);
        %Generate complex noise samples with noise power = 1/SNR_linear
        noise = sqrt(1/(2*SNR_lin)) * (randn(N_R,1) + 1j*randn(N_R,1));
        %Generate a 2x2 MIMO channel
        H = sqrt(1/2)*(randn(N_R,N_T) + 1j*randn(N_R,N_T));
        %Received signal
        y = H*x+noise;
        %Zero-forcing combining--------------------------------------------
        G_ZF = inv(H'*H)*H';
        x_ZF = G_ZF*y;
        bOut_ZF = qamdemod(x_ZF, M,'OutputType','bit','UnitAveragePower',true);
        iterBER_ZF(t)= sum(bOut_ZF ~= b_in)/length(b_in);    
        %MMSE combining----------------------------------------------------
        G_MMSE = inv(H'*H+(1/SNR_lin)*eye(N_R))*H';
        x_MMSE = G_MMSE*y;
        bOut_MMSE = qamdemod(x_MMSE, M,'OutputType','bit','UnitAveragePower',true);
        iterBER_MMSE(t)= sum(bOut_MMSE ~= b_in)/length(b_in);   
        %SVD precoding and combining
        [U,Sigma,V] = svd(H);
        y_precoding = H*V*x+noise;
        x_SVD = U'*y_precoding;
        bOut_SVD = qamdemod(x_SVD, M,'OutputType','bit','UnitAveragePower',true);
        iterBER_SVD(t) = sum(bOut_SVD ~= b_in)/length(b_in);
    end
    avgBER_ZF(i) = mean(iterBER_ZF);
    avgBER_MMSE(i) = mean(iterBER_MMSE);
    avgBER_SVD(i) = mean(iterBER_SVD);
end

figure
semilogy(Eb_N0, avgBER_ZF, 'x-', Eb_N0, avgBER_MMSE, '-o', Eb_N0, avgBER_SVD, '-s')
xlabel('E_b/N_0')
ylabel('BER')
legend('ZF', 'MMSE', 'SVD')
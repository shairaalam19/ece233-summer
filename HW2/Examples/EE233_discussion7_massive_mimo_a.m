clear;

%% System Parameters
r = 2000;               % meters
r_min = 200;            % meters    
rho = 3.4;              % path loss coefficient
K = 20;                 % number of users
Pu_range = 28:-3:10;    % uplink transmit power from UEs in dBm
noise_power = -100;     % noise power at each antenna and UE in dBm
noise_power_lin = 10^(noise_power/10);
X = eye(K);             % uplink pilots for channel estimation
max_iter = 1000;
Pu_lin_range = 10.^(Pu_range/10);

%% Channel estimation
ch_error = zeros(size(Pu_range));
figure(1); hold on;
for M =[100, 500] % number of antennas1
    for iter = 1:max_iter
        % generate channels
        H = zeros(M,K);
        D = zeros(K);
        for k = 1:K
            % generate distance between UE-k and BS uniformly between 200 and 2000
            dk = 200+1800*rand; 
            H(:,k) = sqrt(dk^(-rho)) * sqrt(1/2)*(randn(M,1) + 1j* randn(M,1));
            D(k,k) = dk^(-rho);
        end
        
        for Pu_lin = Pu_lin_range
            % received signal when uplink pilots are transmitted
            noise = sqrt(noise_power_lin/2)* (randn(M,K) + 1j*randn(M,K));
            Y = sqrt(Pu_lin)*H*X + noise;
            
            % LMMSE estimate the channel matrix
            H_est = Y* inv( Pu_lin*D + noise_power_lin * eye(K))* sqrt(Pu_lin)*D;
            
            ch_error(Pu_lin == Pu_lin_range) = ch_error(Pu_lin == Pu_lin_range) + (norm(H-H_est, 'fro'))^2/ (norm(H, 'fro'))^2;
        end
        
    end
    
    ch_error = ch_error/ max_iter;
    plot(Pu_range, ch_error, 'linewidth',3);
    xlabel('Uplink power (dBm)');ylabel('Normalized channel est. error');
    grid on;
end
legend('M=100', 'M=500');
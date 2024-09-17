clear;

%% System Parameters
r = 2000;               % meters
r_min = 200;            % meters    
rho = 3.4;              % path loss coefficient
K = 20;                 % number of users
noise_power = -100;     % noise power at each antenna and UE in dBm
noise_power_lin = 10^(noise_power/10);
max_iter = 1000;
P0 = 40;                % total transmit power in DL from BS (dBm)
P0_lin = 10^(P0/10);
PDk_lin = P0_lin/K;

%% beamforming
% Initialization
M_range = [100, 500];
SINR_ZF_DL = zeros(max_iter, length(M_range));
SINR_MRC_DL = zeros(max_iter, length(M_range));

for M = M_range
    for iter = 1:max_iter
        V_MRC = zeros(M,K); % MRC beamforming vectors;
        V_ZF = zeros(M,K); % ZF beamforming vectors;
        % generate channels
        H = zeros(M,K);
        D = zeros(K);
        for k = 1:K
            dk = 200+1800*rand; % generate distance between UE-k and BS uniformly between 200 and 2000
            H(:,k) = sqrt(dk^(-rho)) * sqrt(1/2)*(randn(M,1) + 1j* randn(M,1));
            D(k,k) = dk^(-rho);
            
            V_MRC(:,k)= conj(H(:,k))/norm(H(:,k));
        end
        
        % ZF beamforming vectors
        G = H*(inv(H'*H));
        for k = 1:K
            V_ZF(:,k) = G(:,k)/ norm(G(:,k));
        end
              
        % SINR for MRC
        DL_int = 0;% interference to UE-1 signal in downlink
        for l = 2:K
            DL_int = DL_int + PDk_lin*abs(H(:,1).'*V_MRC(:,l))^2;
        end
        SINR_MRC_DL(iter, M==M_range) = PDk_lin* norm(H(:,1))^2 /(DL_int + noise_power_lin);
        
        % SINR for zero forcing
        SINR_ZF_DL(iter, M==M_range) = PDk_lin* abs(V_ZF(:,1)'* H(:,1))^2 / noise_power_lin;
    end
end


M = 100;
figure(6); hold on;
[F,X] = ecdf(SINR_ZF_DL(:,M==M_range));
plot(10*log10(X),F);
[F,X] = ecdf(SINR_MRC_DL(:,M==M_range));
plot(10*log10(X),F);
title("M = 100")


M = 500;
figure(7); hold on;
[F,X] = ecdf(SINR_ZF_DL(:,M==M_range));
plot(10*log10(X),F);
[F,X] = ecdf(SINR_MRC_DL(:,M==M_range));
plot(10*log10(X),F);
title("M = 500")


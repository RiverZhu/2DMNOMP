clc; clear; close all;

rng(1);

% signal size parameter
Nx = 64;
My = 32;

NM_num = Nx * My;
% the status of targets
K_targets = 4;
sigma_n = 1;
T = 10;
SNR_dB = 30;   % integrated SNR
gain_allK = sqrt(10 .^ (SNR_dB / 10) * sqrt(sigma_n)) .* exp(1j * 2 * pi * rand(K_targets, T));

omega_x = [0.2;0.4;0.6;0.7]*2*pi;
omega_y = [0.1;0.3;0.5;0.8]*2*pi;
omega_true = [omega_x,omega_y];

ant_idx_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
ant_idx_My = (0 : (My - 1))' - (My - 1) / 2;

array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);
Ztrue = zeros(Nx,My,T);
for t_idx = 1:T
    for tar_idx = 1:K_targets
        Ztrue(:,:,t_idx) = Ztrue(:,:,t_idx)+gain_allK(tar_idx,t_idx)*array_Fun(omega_x(tar_idx,1),Nx)*(array_Fun(omega_true(tar_idx,2),My).');
    end
end
y_matrix = Ztrue + (sigma_n / 2) * (1j * randn(Nx, My,T) + randn(Nx, My,T));

gamma_mnomp = [4, 4];
P_fa = 0.01;
tau = sigma_n * chi2inv((1 - P_fa) ^ (1 / NM_num), 2 * T) / 2;
[omegaList_tau, gainList_tau, y_residue_matrix] =...
MNOMP2D(y_matrix, tau);


figure(1)
subplot(1,2,1)
stem(omega_true(:,1),mean(abs(gain_allK).^2,2),'ro')
hold on
stem(omegaList_tau(:,1),mean(abs(gainList_tau).^2,2),'b+')
legend('true','estimated')
title('frequency in the first dimension')
subplot(1,2,2)
stem(omega_true(:,2),mean(abs(gain_allK).^2,2),'ro')
hold on
stem(omegaList_tau(:,2),mean(abs(gainList_tau).^2,2),'b+')
legend('true','estimated')
title('frequency in the second dimension')






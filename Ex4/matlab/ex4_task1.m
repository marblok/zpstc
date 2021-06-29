function ex4_task1


N_FSD = 13;
tau_N = (N_FSD-1)/2;

Fp1 = 8000; Fp2 = 10000;

M = Fp1 / gcd(Fp1,Fp2) 
L = Fp2 / gcd(Fp1,Fp2)

eps_all = -linspace(-0.5+0.5/L, 0.5-0.5/L, L);

figure(1)
ha1 = subplot(2,1,1);
h_all = zeros(1, N_FSD*L);
n_all = 0:length(h_all)-1; n_all = n_all - (length(n_all)-1)/2;
for ind = 1:L,
  epsilon = eps_all(ind);
  h = FSDfilter('Lagr', N_FSD, tau_N+epsilon);
%   h = FSDfilter('optimal2', N_FSD, tau_N+epsilon, 0.4);
  h_all(ind:L:end) = h;
  set(stem(ha1, n_all, h_all), 'Marker', 'none');
  pause(0.1)
  drawnow
end

K = 16*8192;
F = linspace(0, L*Fp1, K);

H_all = abs(fft(h_all/L, K));
figure(11)
subplot(2,1,1)
plot(h_all);
subplot(2,1,2)
plot(F, 20*log10(H_all));
set(gca, 'Xlim', [0, L*Fp1/2], 'Ylim', [-100, 3])


dane.h{1}  = h_all;
dane.h{2}  = [L, M];
dane.Fp = Fp1;
save_filter_coef('ex4_task2_h_FSD_all', dane, 1);


function test

Fsymb = 3000;
K = 32;
[x, Fp] = readaudiofile('../ex6_task1.flt', 'cplx');
[xb, Fp] = readaudiofile('../ex6_task1b.flt', 'cplx');

length(x)
Fp

% y = [];
% while length(xb) >= K
%   y = [y; xb(K:-1:1)];
%   xb(1:K) = [];
% end
% xb = [y; xb];

figure(1)
subplot(2,1,1)
plot(real(x)) 
hold on
plot(real(xb), 'r') 
hold off
subplot(2,1,2)
plot(imag(x)) 

K_psd = 8*1024; 

figure(2)
specgram(x)

figure(3)
% psd(x)
F = linspace(-Fp/2, Fp/2, 2*K_psd+1); F(end) = [];
% pwelch(x, blackman(K_psd), ceil(0.5*K_psd), F, Fp);
[Pxx, F] = pwelch(x, blackman(K_psd), ceil(0.5*K_psd), F, Fp);
[Pxx_b, F] = pwelch(xb, blackman(K_psd), ceil(0.5*K_psd), F, Fp);
Pxx_dB = 10*log10(Pxx);
Pxx_b_dB = 10*log10(Pxx_b);
Pxx_dB_max = ceil(max(Pxx_dB)/3)*3; Y_range = [-100, 0]+Pxx_dB_max;
plot([0,0], Y_range, 'r--');
set(gca, 'Xlim', [-Fp/2, Fp/2]/1000, 'Ylim', Y_range);
hold on
for ind = 1:ceil(K/2),
  plot(ind*Fsymb/1000*[1,1], Y_range, 'r--');
  plot(-ind*Fsymb/1000*[1,1], Y_range, 'r--');

  plot((ind-0.5)*Fsymb/1000*[1,1], Y_range, 'k-');
  plot(-(ind-0.5)*Fsymb/1000*[1,1], Y_range, 'k-');
end
plot(F/1000,Pxx_dB);
plot(F/1000,Pxx_b_dB, 'm:');
hold off
% set(gca,'Xtick', (-(floor(K/2)+1/2):1/2:(ceil(K/2)+1/2))*Fsymb/1000, 'Xgrid', 'on');

xlabel('Frequency [kHz]')
ylabel('Power spectral density [dB/Hz]')
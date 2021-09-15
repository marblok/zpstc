function test_FMD_signal(filename, N_DFT, Fsymb, use_only_real_part)

if nargin == 0,
  % parametry sygna³u wielotonowego
  Fsymb = 3000;
  N_DFT = 16;
  filename = '../ex6_task1.flt';
  use_only_real_part = 1;
end

try
  [x, Fp] = readaudiofile(filename, 'cplx');
catch
  [x, Fp] = readaudiofile(filename);
end

K_psd = 8*1024; 
K_spect = 1024; 

if use_only_real_part ~= 0,
  x = real(x);
  % freq_range = 'onesided';
  if rem(K_spect,2) == 1,
    K_spect = K_spect+1;
  end
  F = linspace(0, Fp/2, K_spect+1); F(end) = [];
else
  % freq_range = 'centered';
  F = linspace(-Fp/2, Fp/2, 2*K_spect+1); F(end) = [];
end

length(x)
Fp

figure(1)
subplot(2,1,1)
plot(real(x)) 
subplot(2,1,2)
plot(imag(x)) 

F_psd = linspace(-Fp/2, Fp/2, 2*K_psd+1); F_psd(end) = [];

figure(12)
set(12, 'Name', 'Spectrogram analizowanego sygna³u');
% specgram(x)
spectrogram(x, blackman(K_spect), ceil(0.5*K_spect), F, Fp, 'yaxis');
colormap(jet)

figure(13)
set(13, 'Name', 'Widmo gêstoœci mocy analizowanego sygna³u');
% psd(x)
% pwelch(x, blackman(K_psd), ceil(0.5*K_psd), F, Fp);
Pxx = pwelch(x, blackman(K_psd), ceil(0.5*K_psd), F_psd, Fp);
Pxx_dB = 10*log10(Pxx);
Pxx_dB_max = ceil(max(Pxx_dB)/3)*3; Y_range = [-100, 0]+Pxx_dB_max;
plot([0,0], Y_range, 'r--');
set(gca, 'Xlim', [-Fp/2, Fp/2]/1000, 'Ylim', Y_range);
hold on
for ind = 1:ceil(N_DFT/2),
  plot(ind*Fsymb/1000*[1,1], Y_range, 'r--');
  plot(-ind*Fsymb/1000*[1,1], Y_range, 'r--');

  plot((ind-0.5)*Fsymb/1000*[1,1], Y_range, 'k-');
  plot(-(ind-0.5)*Fsymb/1000*[1,1], Y_range, 'k-');
end
plot(F_psd/1000,Pxx_dB);
hold off
% set(gca,'Xtick', (-(floor(K/2)+1/2):1/2:(ceil(K/2)+1/2))*Fsymb/1000, 'Xgrid', 'on');

xlabel('Frequency [kHz]')
ylabel('Power spectral density [dB/Hz]')


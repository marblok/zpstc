function test_task1_signal(mode,matched, offset_0, phase_offset)
% test_task1_signal('zad1', 'srRC', 0, -2*pi/32)
% test_task1_signal('zad1b', 'srRC', 0, -2*pi/32)
%
% offset_0 - przesuni�cie w pr�bkach chwili pr�bkowania sybmoli
% phase_offset - korekta fazy pocz�tkowej heterodynny

channel_no = 1;
% channel_no = 4;
% channel_no = 6;


% if nargin == 0,
%   ind = 2*(3+2*4);
% end

if nargin == 0
  % matched = 'srRC' 
  matched = 'rect'
  % mode = 'zad1'
  mode = 'zad1'
  
  offset_0 = 0;
  phase_offset = 0;
elseif  nargin == 2
  offset_0 = 0;
  phase_offset = 0;
elseif  nargin == 3
  phase_offset = 0;
end



N_DFT = 32; M = 8;
N_symb = N_DFT;

switch mode
  case 'zad1',
    [x, Fp] = fileread('../ex6_task1.flt', 'cplx');
    offset = 5*N_symb + offset_0;
  case 'zad1b',
    [x, Fp] = fileread('../ex6_task1b.flt', 'cplx');
    offset = 6*N_symb + offset_0;
end
% % [y_ref, Fp] = fileread('../ex6_task1a.flt', 'cplx');

[x_ch1, F_symb] = fileread('../ex6_task1_ch1.flt', 'cplx');


switch matched,
  case 'rect',
    h_matched = ones(1,N_DFT); % du�y poziom interferencji ISI/ICI
    offset = offset + 0.5*N_symb;

  case 'srRC',
    [coef, ver] = load_filter_coef('../matlab/ex6_task1_h_rc.coef')
    h_rc = coef.h;
    F_symb = coef.Fp;

    h_matched = N_symb*h_rc;
    offset = offset + 5*N_symb;
end
N_matched = length(h_matched);

K_psd = 1024;
N_symb = N_DFT

figure(11)
subplot(2,1,1)
n = 0:N_matched-1;
plot(n-(N_matched-1)/2, h_matched)
hold on
n = 0:(2*N_matched-1)-1;
plot(n-2*(N_matched-1)/2, conv(h_matched,h_matched), 'r')
plot([-N_matched, N_matched], [0, 0], 'k')
hold off
set(gca,'xlim', [-N_matched, N_matched], 'xtick', N_symb*[-16:16], 'xgrid', 'on')

subplot(2,1,2)
[H, F_h] = freqz(h_matched, 1, K_psd, F_symb*N_DFT);
[H2, F_h] = freqz(conv(h_matched, h_matched), 1, K_psd, F_symb*N_DFT);
plot(F_h, 20*log10(abs(H)))
hold on
plot(F_h, 20*log10(abs(H2)), 'r')
hold off
set(gca, 'xlim', [0, F_symb/2*N_DFT], 'xtick', [0:N_DFT]*F_symb/2, 'xgrid', 'on')
set(gca, 'Ylim', [-170, 3])

% channel_no = 1;

% size(x)
% size(y)
% figure
% plot(x-y)
% return

Fsymb = Fp/N_DFT
% K = 32; M = 8;

figure(1)
subplot(4,1,1)
plot(real(x), 'b')
title('wygenerowany sygna� wielokana�owy')
%axis equal
subplot(4,1,2)
plot(imag(x), 'r')
%axis equal

subplot(2,1,2)
F_psd = linspace(-Fp/2, Fp/2, 2*K_psd+1); F_psd(end) = [];
% psd(het.*y)
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
k = [-ceil(N_DFT/2):ceil(N_DFT/2)]';
ticks = k*Fsymb/1000;
set(gca, 'xtick', ticks, 'xticklabel', num2str(k))


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
n = [0:length(x)-1].';
het = exp(-j*channel_no*(2*pi)/N_DFT*n + j*phase_offset);
y = het .* x;

figure(2)
subplot(4,1,1)
plot(real(y), 'b')
title('sygna� heterodynowany w d�')
%axis equal
subplot(4,1,2)
plot(imag(y), 'r')
%axis equal

subplot(2,1,2)
Pxx = pwelch(y, blackman(K_psd), ceil(0.5*K_psd), F_psd, Fp);
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
k = [-ceil(N_DFT/2):ceil(N_DFT/2)]';
ticks = k*Fsymb/1000;
set(gca, 'xtick', ticks, 'xticklabel', num2str(k))

% y = y(21+1:end);
y2 = filter(h_matched, 1, y);

figure(3)
subplot(4,1,1)
plot(real(y2), 'b')
title('sygna� heterodynowany w d� po filtrze dopasowanym')
%axis equal
subplot(4,1,2)
plot(imag(y2), 'r')
%axis equal

subplot(2,1,2)
Pxx = pwelch(y2, blackman(K_psd), ceil(0.5*K_psd), F_psd, Fp);
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
k = (-ceil(N_DFT/2):ceil(N_DFT/2))';
ticks = k*Fsymb/1000;
set(gca, 'xtick', ticks, 'xticklabel', num2str(k))

m = 0:length(x_ch1)-1;
figure(3)
subplot(4,1,1)
hold on
plot(N_DFT*m+offset, real(x_ch1), 'k');
hold off
subplot(4,1,2)
hold on
plot(N_DFT*m+offset, imag(x_ch1), 'k');
hold off

symbols = y2(offset:N_symb:end);
figure(4)
plot(symbols, 'b.')
axis equal

return

% figure(2)
% subplot(2,1,1)
% stem(real(x)/max(abs(x)), 'b')
% hold on
% stem(real(y_ref), 'r')
% hold off
% subplot(2,1,2)
% [c, l] = xcorr(real(y), real(y_ref));
% plot(l, c)
% 
% figure(3)
% subplot(2,1,1)
% stem(imag(y)/max(abs(y)), 'b')
% hold on
% stem(imag(y_ref), 'r')
% hold off
% subplot(2,1,2)
% [c, l] = xcorr(imag(y), imag(y_ref));
% plot(l, c)
% 
% % het = 1;
% figure(13)
% 
% % % return
% % xt{1} = ((-1).^(1:length(y)).').*y/max(abs(y));
% xt{1} = het.*y/max(abs(y));
% xt{2} = y_ref;
% trajekt_(xt, 10, 1:M:length(y), 0, 1)
% tmp = xt{1};
% eyediagram2(tmp(1:100*M), 0, M, 15);
% pause
% 
% y = y(1:M:end);
% % y = y.*exp(j*pi/4);
% s_re = real(y) > 0;
% s_im = imag(y) > 0;
% 
% bin = zeros(1,2*length(s_re));
% bin(1:2:end) = s_re;
% bin(2:2:end) = s_im;
% 
% bits_per_char = 8;
% %ind = 0;
% 
% h = fopen('../Ex6_task1.cpp');
% bin2 = double(fread(h, 800, 'ubit1')).';
% fclose(h)
% % 
% figure(4)
% subplot(2,1,1)
% [c, l] = xcorr(s_re, bin2(1:2:end));
% plot(l, c)
% subplot(2,1,2)
% [c, l] = xcorr(s_re, bin2(2:2:end));
% plot(l, c)
% 
% % return  
% figure(5)
% subplot(2,1,1)
% stem(bin, 'b')
% hold on
% stem(bin2, 'r')
% hold off
% subplot(2,1,2)
% [c, l] = xcorr(bin, bin2);
% plot(l, c)
% 
% % bin = bin_root.'
% bin = bin((ind+1):end);
% % for ind= 1:bits_per_char,
% for ind= 1,
%   bin_ = bin(ind:end);
%   B = floor(length(bin_)/bits_per_char)*bits_per_char;
%   
%   tekst = bin_(1:bits_per_char:B);
%   for bit = 1:bits_per_char-1,
% %     tekst = tekst + (2^bit)*bin_((bit+1):bits_per_char:B);
%     tekst = tekst + (2^((bits_per_char-1)-bit))*bin_((bit+1):bits_per_char:B);
%   end
%   char(tekst)
% end
% 
% 
% % 
% % B_pre = length(pre)*bits_per_char;
% % bin_pre = zeros(B_pre,1);
% % for ind = 0:bits_per_char-1,
% %   bin_pre((ind+1):bits_per_char:end) = rem(floor(pre / (2^ind)),2); 
% % end

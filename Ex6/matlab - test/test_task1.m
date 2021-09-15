function test_task1(ind)
if nargin == 0,
  ind = 2*(3+2*4);
end
[x, Fp] = fileread('../ex6_task1.flt', 'cplx');
[y, Fp_y] = fileread('../ex6_task1b.flt', 'cplx');
[y_ref, Fp] = fileread('../ex6_task1a.flt', 'cplx');

channel_no = 1;

% size(x)
% size(y)
% figure
% plot(x-y)
% return

N_DFT = 32; M = 8;
Fsymb = Fp_y/N_DFT
% K = 32; M = 8;

figure(1)
subplot(2,1,1)
plot(x)
axis equal
subplot(2,2,3)
plot(y, 'b')
axis equal
subplot(2,2,4)
plot(y_ref, 'b.')
axis equal
% % plot(y)

% y = y(21+1:end);

figure(2)
subplot(2,1,1)
stem(real(y)/max(abs(y)), 'b')
hold on
stem(real(y_ref), 'r')
hold off
subplot(2,1,2)
[c, l] = xcorr(real(y), real(y_ref));
plot(l, c)

figure(3)
subplot(2,1,1)
stem(imag(y)/max(abs(y)), 'b')
hold on
stem(imag(y_ref), 'r')
hold off
subplot(2,1,2)
[c, l] = xcorr(imag(y), imag(y_ref));
plot(l, c)

n = [0:length(y)-1].';
het = exp(-j*channel_no*(2*pi)/N_DFT*n);
% het = 1;
figure(13)
K_psd = 1024;
F_psd = linspace(-Fp_y/2, Fp_y/2, 2*K_psd+1); F_psd(end) = [];
% psd(het.*y)
Pxx = pwelch(het.*y, blackman(K_psd), ceil(0.5*K_psd), F_psd, Fp_y);
Pxx_dB = 10*log10(Pxx);
Pxx_dB_max = ceil(max(Pxx_dB)/3)*3; Y_range = [-100, 0]+Pxx_dB_max;
plot([0,0], Y_range, 'r--');
set(gca, 'Xlim', [-Fp_y/2, Fp_y/2]/1000, 'Ylim', Y_range);
hold on
for ind = 1:ceil(N_DFT/2),
  plot(ind*Fsymb/1000*[1,1], Y_range, 'r--');
  plot(-ind*Fsymb/1000*[1,1], Y_range, 'r--');

  plot((ind-0.5)*Fsymb/1000*[1,1], Y_range, 'k-');
  plot(-(ind-0.5)*Fsymb/1000*[1,1], Y_range, 'k-');
end
plot(F_psd/1000,Pxx_dB);
hold off

% % return
% xt{1} = ((-1).^(1:length(y)).').*y/max(abs(y));
xt{1} = het.*y/max(abs(y));
xt{2} = y_ref;
trajekt_(xt, 10, 1:M:length(y), 0, 1)
tmp = xt{1};
eyediagram2(tmp(1:100*M), 0, M, 15);
pause

y = y(1:M:end);
% y = y.*exp(j*pi/4);
s_re = real(y) > 0;
s_im = imag(y) > 0;

bin = zeros(1,2*length(s_re));
bin(1:2:end) = s_re;
bin(2:2:end) = s_im;

bits_per_char = 8;
%ind = 0;

h = fopen('../Ex6_task1.cpp');
bin2 = double(fread(h, 800, 'ubit1')).';
fclose(h)
% 
figure(4)
subplot(2,1,1)
[c, l] = xcorr(s_re, bin2(1:2:end));
plot(l, c)
subplot(2,1,2)
[c, l] = xcorr(s_re, bin2(2:2:end));
plot(l, c)

% return  
figure(5)
subplot(2,1,1)
stem(bin, 'b')
hold on
stem(bin2, 'r')
hold off
subplot(2,1,2)
[c, l] = xcorr(bin, bin2);
plot(l, c)

% bin = bin_root.'
bin = bin((ind+1):end);
% for ind= 1:bits_per_char,
for ind= 1,
  bin_ = bin(ind:end);
  B = floor(length(bin_)/bits_per_char)*bits_per_char;
  
  tekst = bin_(1:bits_per_char:B);
  for bit = 1:bits_per_char-1,
%     tekst = tekst + (2^bit)*bin_((bit+1):bits_per_char:B);
    tekst = tekst + (2^((bits_per_char-1)-bit))*bin_((bit+1):bits_per_char:B);
  end
  char(tekst)
end


% 
% B_pre = length(pre)*bits_per_char;
% bin_pre = zeros(B_pre,1);
% for ind = 0:bits_per_char-1,
%   bin_pre((ind+1):bits_per_char:end) = rem(floor(pre / (2^ind)),2); 
% end

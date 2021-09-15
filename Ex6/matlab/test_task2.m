function test_task2(ind)
if nargin == 0,
  ind = 2*(3+2*4);
end
[x, Fp] = readaudiofile('../ex6_task1.flt');
[y, Fp] = readaudiofile('../ex6_task2a.flt', 'cplx');
[y_ref, Fp] = readaudiofile('../ex6_task1a.flt', 'cplx');

K = 32; M = 8;
% K = 16; M = 8;

figure(1)
subplot(2,1,1)
plot(x)
subplot(2,2,3)
plot(y, 'b.')
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
het = exp(-j*(8)*pi/(32)*n);
% het = 1;
figure(13)
psd(het.*y)
% return

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
bin2 = fread(h, 800, 'ubit1').';
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

function test_task2(symb_offset,phase_offset, bit_offset)
% test_task2(2+4*ile, 0) 
% test_task2(10, 0) - % pe�na kompensacja op�nienia
% demoulacja pierwszego podkana�u
% test_task2(10,0,2) dla ex6_task1.flt
% test_task2(11,0,2) dla ex6_task1b.flt
if nargin == 0,
%   ind = 2*(3+2*4);
  symb_offset = 10;
  phase_offset = 0;
  bit_offset = 2;
end
phasor = exp(j*phase_offset);

% dla ex6_task1.flt
[x, Fp] = fileread('../ex6_task1.flt');
[y, Fp] = fileread('../ex6_task2_symb_A_a.flt', 'cplx');
% dla ex6_task1b.flt
[x, Fp] = fileread('../ex6_task1b.flt');
[y, Fp] = fileread('../ex6_task2_symb_B_a.flt', 'cplx');

[y_ref, Fp] = fileread('../ex6_task1a.flt', 'cplx');

y = y(symb_offset+1:end);

% y = y_ref; % - przesuni�cie o 4


mode = 'QPSK_A', % 00 - 1, 10 - j, 01 - -j, 11 - -1

% % y = y  - (0.5+j*0.5);
% y_ref = (real(y_ref) > 0) + j*(imag(y_ref) > 0);
K = 32; M = 32;
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

% return 

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
% het = exp(-j*(8)*pi/(32)*n);
het = 1;
% het = 1;
figure(13)
psd(het.*y)
% return

if M < K
  % xt{1} = ((-1).^(1:length(y)).').*y/max(abs(y));
  xt{1} = het.*y/max(abs(y));
  xt{2} = y_ref;
  trajekt_(xt, 10, 1:M:length(y), 0, 1)
  tmp = xt{1};
  eyediagram2(tmp(1:100*M), 0, M, 15);
  pause

  y = y(1:M:end);
end
y = y * phasor;


switch mode,
  case 'QPSK_A', % �le ? % 00 - 1, 10 - j, 01 - -j, 11 - -1
    % 00 - 1, 01 - j, 10 - -j, 11 - -1
    for n = 1:length(y),
      if real(y(n)) > abs(imag(y(n))), % 1
        s_re(n) = 0.0;
        s_im(n) = 0.0;
      elseif -real(y(n)) > abs(imag(y(n))), % -1
        s_re(n) = 1.0;
        s_im(n) = 1.0;
      elseif imag(y(n)) > abs(real(y(n))),  % j
        s_re(n) = 1.0;   s_im(n) = 0.0;
%         s_re(n) = 0.0;   s_im(n) = 1.0;
      elseif -imag(y(n)) > abs(real(y(n))), % -j
        s_re(n) = 0.0; s_im(n) = 1.0;
%         s_re(n) = 1.0; s_im(n) = 0.0;
      else
        s_re(n) = 0.0;
        s_im(n) = 0.0;
      end
    end
    
  case 'QPSK_B', % ???
    s_re = real(y) > 0;
    s_im = imag(y) > 0;
end


bin = zeros(1,2*length(s_re));
bin(1:2:end) = s_re;
bin(2:2:end) = s_im;

bits_per_char = 8;
%ind = 0;

h = fopen('../Ex6_task1.cpp');
bin2 = fread(h, 800, 'ubit1', 'ieee-be').'; % ???
fclose(h);
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
stem(bin2, 'rx')
hold off
subplot(2,1,2)
[c, l] = xcorr(bin, bin2);
plot(l, c)

% % bin = bin_root.'
% bin = bin((ind+1):end);
% for ind= 1:bits_per_char,
for ind= bit_offset+1,
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

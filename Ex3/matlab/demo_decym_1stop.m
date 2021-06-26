function demo_decym_1stop
% Decymacja jednostopniowa
close all;
clear all;
 
%definiowanie wejœciowej i wyjœciowej szybkoœci próbkowania
fold = 10240;
fnew = 320;
 
%obliczanie krotnoœci decymacji
M = fold/fnew;
 
%obliczanie granic pasma przepustowego i pasma zaporowego
vp = 0.4/M;
vz = 0.5/M;
vg = 0.5/M;
vp1 = 0.4/M;
vz1 = 0.6/M;
 
Ap = 0.2;
dp = -((1-10^(Ap/20))/(1+10^(Ap/20)));
As = 60;
ds = 10^(-As/20);
 
%definiowanie sygna³u testowego
n = 0:5*fold-1;
x = sin(2*pi*9/fold*n)+1/4*sin((2*pi*27/fold*n)+pi/4)+1/16*sin((2*pi*36/fold*n)+pi/6)+sin(2*pi*49/fold*n);
x = x + randn(size(x))/100000; 
 
%Widmo sygna³u wejsciowego
figure(1)
plot(n,20*log10(abs(fft(x))));
axis ([0 fold/2 -100 100]);
 
%projektowanie filtru
[nn,fo,mo,w] = remezord( [vp vz], [1 0], [dp ds], 1);
 
h = remez(nn,fo,mo,w);
[H, F2] = freqz(h, 1, 8192*32, 1);
 
%Charakterystyka amplitudowa filtru H'
figure(2)
subplot(2,2,1);
plot(F2,20*log10(abs(H)));
hold on;
plot ([vp vp],[-80 5],'k--');
plot ([vz vz],[-80 5],'k:');
plot ([vg vg],[-80 5],'k-');
axis ([0 0.5 -80 5]);
title('a');
 
 
y = conv(h,x);
 
p = length(h)-1;
yyy=y(p+1:M:end-1);
 
wind = blackman(length(yyy));
Y = 20*log10(abs(fft(wind'.*yyy,32768)));
 
%Widmo sygnalu x[n] po filracji
subplot(2,2,2);
plot(Y);
title('b');
axis ([0 length(Y)/2 -100 55]);
 
%projektowanie filtru2
[nn2,fo2,mo2,w2] = remezord( [vp1 vz1], [1 0], [dp ds], 1);
 
h2 = remez(nn2,fo2,mo2,w2);
[H2, F22] = freqz(h2, 1, 8192*32, 1);
 
%Charakterystyka amplitudowa filtru H2
subplot(2,2,3);
plot(F22,20*log10(abs(H2)));
hold on;
plot ([vp1 vp1],[-80 5],'k--');
plot ([vz1 vz1],[-80 5],'k:');
plot ([vg vg],[-80 5],'k-');
axis ([0 0.5 -80 5]);
title('c');
 
y2 = conv(h2,x);
 
p2 = length(h2)-1;
yyy2=y2(p2+1:M:end-1);
 
n22 = 0:fold-1;
wind2 = blackman(length(yyy2));
Y2 = 20*log10(abs(fft(wind2'.*yyy2,32768)));
 
%Widmo sygnalu x[n] po filtracji
subplot(2,2,4);
plot(Y2);
title('d');	
axis ([0 length(Y2)/2 -100 55]);
 
 
 
K = 8192*32; % zeropadding
wind11 = blackman(length(yyy));
Y11 = 20*log10(abs(fft(wind11'.*yyy,K)));
wind22 = blackman(length(yyy2));
Y22 = 20*log10(abs(fft(wind22'.*yyy2,K)));
 
%Widmo sygnalu x[n] po decymacji
figure(3)
subplot(2,1,1);
 
plot(Y11);
title('a');
axis ([0 length(Y11)/2 -100 50]);
subplot(2,1,2);
plot(Y22);
title('b');
axis ([0 length(Y22)/2 -100 50]);
 
for i=1:3,
    set(i, 'color', 'w');
end

function demo_inter_1stop
close all;
clear all;
 
%definiowanie wejœciowej i wyjœciowej szybkoœci próbkowania
fold = 100;
fnew = 400;
 
%obliczanie krotnoœci interpolacji
L = fnew/fold;
 
%obliczanie granic pasma przepustowego i pasma zaporowego
vp = 0.4/L;
vz = 0.5/L;
vg = 0.5/L;
vp1 = 0.4/L;
vz1 = 0.6/L;
 
Ap = 0.2;
dp = -((1-10^(Ap/20))/(1+10^(Ap/20)));
As = 60;
ds = 10^(-As/20);
 
%definiowanie sygna³u testowego
n = 0:5*fold-1;
x = sin(2*pi*9/fold*n)+1/4*sin((2*pi*27/fold*n)+pi/4)+1/16*sin((2*pi*36/fold*n)+pi/6)+sin(2*pi*49/fold*n);
x = x + randn(size(x))/100000; 
 
figure(1)
subplot(1,2,1);
plot(n,20*log10(abs(fft(x))));
title('a');
axis ([0 length(x)/2 -100 50]);
 
%uzupe³nianie zerami
xx = zeros(1,L*length(x));
xx(1:L:length(x)*L)=x;
 
n2 = 0:fnew-1;
 
 
subplot(1,2,2);
plot(n2,20*log10(abs(fft(xx(1:L*fold)))));
title('b');
axis ([0 length(x)/2 -100 50]);
 
%projektowanie filtru
[nn,fo,mo,w] = remezord( [vp vz], [1 0], [dp ds], 1);
h = L*remez(nn,fo,mo,w);
[H, F2] = freqz(h/L, 1, 8192*32, 1);
 
 
figure(2)
subplot(2,2,1);
plot(F2,20*log10(abs(H)));
hold on;
plot ([vp vp],[-80 5],'k--');
plot ([vz vz],[-80 5],'k:');
plot ([vg vg],[-80 5],'k-');
axis ([0 0.5 -80 5]);
title('a');
 
 
y = conv(h,xx);
 
p = length(h)-1;
yy=y(p+1:1:end);
 
 
 
wind = blackman(length(yy));
Y = 20*log10(abs(fft(wind'.*yy,32768)));
 
subplot(2,2,2);
plot(Y);
axis ([0 length(Y)/2 -155 50]);
title('b');
 
%projektowanie filtru2
[nn2,fo2,mo2,w2] = remezord( [vp1 vz1], [1 0], [dp ds], 1);
h2 = L*remez(nn2,fo2,mo2,w2);
[H2, F22] = freqz(h2/L, 1, 8192*32, 1);
 
subplot(2,2,3);
plot(F22,20*log10(abs(H2)));
hold on;
plot ([vp1 vp1],[-80 5],'k--');
plot ([vz1 vz1],[-80 5],'k:');
plot ([vg vg],[-80 5],'k-');
axis ([0 0.5 -80 5]);
title('c');
 
 
y2 = conv(h2,xx);
 
p2 = length(h2)-1;
yy2=y2(p2+1:1:end);
 
 
wind2 = blackman(length(yy2));
Y2 = 20*log10(abs(fft(wind2'.*yy2,32768)));
 
subplot(2,2,4);
plot(Y2);
%plot(20*log10(abs(fft(yy2))));
axis ([0 length(Y2)/2 -155 50]);
title('d');
 
 
figure(3);
K = 8192*32; % zeropadding
okno = chebwin(length(yy), 150).';
 
YY= 20*log10(abs(fft(yy.*okno, K)));
plot(YY);
 
hold on
okno = chebwin(length(yy2), 150).';
 
plot(20*log10(abs(fft(yy2.*okno, K))), 'r');
axis ([0 length(YY) -100 50]);
hold off
 
for i=1:3,
    set(i, 'color', 'w');
end


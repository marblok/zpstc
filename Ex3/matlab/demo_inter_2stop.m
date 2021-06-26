function demo_inter_2stop
close all;
clear all;
 
%definiowanie wejœciowej i wyjœciowej szybkoœci próbkowania
fold = 320;
fnew = 10240;
Fg = 128;
 
%obliczanie krotnoœci interpolacji
L = fnew/fold;

% implementacja jednostopniowa
vp = Fg/fnew
vz = (fold/2)/fnew
 
Ap = 0.1;
dp = -((1-10^(Ap/20))/(1+10^(Ap/20)));
As = 60;
ds = 10^(-As/20);

%projektowanie filtru
[nn,fo,mo,w] = remezord( [vp vz], [1 0], [dp ds], 1);
h = L*remez(nn,fo,mo,w);
[H, F] = freqz(h/L, 1, 8192*32, fnew);

figure(1)
subplot(2,1,1);
plot(F,20*log10(abs(H)));
hold on;
plot ([vp vp]*fnew,[-80 5],'k--');
plot ([vz vz]*fnew,[-80 5],'k:');
hold off;
axis ([0 0.5*fnew -80 5]);

N = length(h)


L1 = 8;
%obliczanie granic pasma przepustowego i pasma zaporowego
% filtr w pierwszym stopniu pracuje z szybkoœci¹ L1*fold)
vp = Fg/(L1*fold)
vz = (fold/2)/(L1*fold)
% % vz = 0.5/L;
% vg = 0.5/L;
% vp1 = 0.4/L;
% vz1 = 0.6/L;
 
Ap = 0.05;
dp = -((1-10^(Ap/20))/(1+10^(Ap/20)));
As = 60;
ds = 10^(-As/20);
 
 
%projektowanie filtru
[nn,fo,mo,w] = remezord( [vp vz], [1 0], [dp ds], 1);
h1 = L1*remez(nn,fo,mo,w);
[H1, F2] = freqz(h1/L1, 1, 8192*32, 1);
 
N1 = length(h1)
 
figure(2)
subplot(2,2,1);
plot(F2,20*log10(abs(H1)));
hold on;
plot ([vp vp],[-80 5],'k--');
plot ([vz vz],[-80 5],'k:');
hold off;
axis ([0 0.5 -80 5]);
title('a');

subplot(2,2,2);
plot(F2,20*log10(abs(H1)));
hold on;
plot ([vp vp],[-80 5],'k--');
plot ([vz vz],[-80 5],'k:');
hold off;
axis ([0 vz -Ap Ap]);
title('b');
 

% drugi stopieñ interpolacji
L2 = 4;
% filtr pracuje z szybkoœci¹ fnew
vp1 = Fg/fnew;
vz1 = (L1*fold-Fg)/fnew

%projektowanie filtru2
[nn2,fo2,mo2,w2] = remezord( [vp1 vz1], [1 0], [dp ds], 1);
h2 = L2*remez(nn2,fo2,mo2,w2);
[H2, F22] = freqz(h2/L2, 1, 8192*32, 1);
 
subplot(2,2,3);
plot(F22,20*log10(abs(H2)));
hold on;
plot ([vp1 vp1],[-80 5],'k--');
plot ([vz1 vz1],[-80 5],'k:');
% plot ([vg vg],[-80 5],'k-');
hold off;
axis ([0 0.5 -80 5]);
title('c');
 

subplot(2,2,4);
plot(F2,20*log10(abs(H2)));
hold on;
plot ([vp1 vp1],[-80 5],'k--');
plot ([vz1 vz1],[-80 5],'k:');
hold off;
axis ([0 vz/L2 -Ap Ap]);
title('d');
 

N2 = length(h2)


h1_L2 = zeros(N1*L2-(L2-1), 1);
h1_L2(1:L2:end) = h1;
[H1_L2, F1_L2] = freqz(h1_L2/L1, 1, 8192*32, 1);

h_all = conv(h1_L2, h2);
N_all = length(h_all)

figure(3)
subplot(2,1,1)
plot(F1_L2,20*log10(abs(H1_L2)));
hold on;
plot(F2,20*log10(abs(H2)), 'r');
hold off
axis ([0 0.5 -80 5]);
title('charakterystyki filtrów sk³adowych przeniesionych na szybkoœæ próbkowania fnew');

subplot(2,1,2)
plot(F1_L2,20*log10(abs(H1_L2.*H2)));

axis ([0 0.5 -80 5]);
title('charakterystyka zbiorcza');

figure(1)
subplot(2,1,2)
plot(F1_L2*fnew,20*log10(abs(H1_L2.*H2)));
axis ([0 0.5*fnew -80 5]);

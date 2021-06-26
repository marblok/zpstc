function N = id_1(nr_sygnalu, fold, fnew, Li, Mi, FPi, FSi)

L = fnew / gcd(fold, fnew);
M = fold / gcd(fold, fnew);

L1 = Li(1);
M1 = Mi(1);

fp1 = FPi(1);
fs1 = FSi(1);


% %definiowanie wejœciowej i wyjœciowej szybkoœci próbkowania
% fold = 320;
% fnew = 10240;
% 
% %obliczanie krotnoœci
% L = fnew/gcd(fold,fnew);
% M = fold/gcd(fold,fnew);
% L1 = L;
% M1 = M;

%obliczanie granic pasma przepustowego i pasma zaporowego
% fg1 = min(0.5/L,0.5/M);

%definiowanie sygna³u testowego
x = gen_signal(nr_sygnalu, fold, L, M);

%obliczanie zafalowan w pasmach przepustowym i zaporowym
Ap = 0.1;
dp = -((1-10^(Ap/20))/(1+10^(Ap/20)));
As = 60;
ds = 10^(-As/20);

%uzupe³nianie zerami(interpolacja – wstawiamy L1-1 zer pomiedzy ka¿d¹ parê próbek)
xx = zeros(1,L1*length(x));
xx(1:L1:length(x)*L1)=x;

%projektowanie filtru
[nn,fo,mo,w] = firpmord( [fp1 fs1], [1 0], [dp ds], 1);
if nn < 3, nn=3; end;

if nn > 1000,
  N = nn;
  disp('wymagany rz¹d filtru > 1000')
  return
end
tic
h1 = L1*firpm(nn,fo,mo,w);
toc
[H, F2] = freqz(h1/L1, 1, 8192*32, fold*L);


%splot sygna³u po interpolacji z odpowiedzi¹ impulsow¹ filtru
y = conv(h1,xx);

% usuwanie stanu przejsciowego i decymacja sygba³u - wybieramy co M-ta
% probke
p = length(h1)-1;
y11 = y(p+1:M1:end);



%wykresy

%widmo sygna³u wejœciowego
wind0 = blackman(length(x));
Y0 = 20*log10(abs(fftshift(fft(wind0'.*x,32768))));
F_Y0 = linspace(-fold/2, fold/2, 32768+1); F_Y0(end) = [];
figure(2)
subplot(1,2,1);
plot(F_Y0, Y0);
set(gca, 'xlim', [-fold/2, fold/2]);
xlabel('F [Hz]');
title('syg. wejsciowy');

%widmo sygba³u wyjsciowego
figure(2)
subplot(1,2,2);
wind1 = blackman(length(y11));
Y1 = 20*log10(abs(fftshift(fft(wind1'.*y11,32768))));
F_Y1 = linspace(-fnew/2, fnew/2, 32768+1); F_Y1(end) = [];
plot(F_Y1, Y1);
set(gca, 'xlim', [-fnew/2, fnew/2]);
xlabel('F [Hz]');
title('syg. wyjsciowy');



%charakterystyka amplitudowa filtru 1
figure(3)
plot(F2,20*log10(abs(H)));
hold on;
plot ([fp1 fp1]*fold*L,[-80 5],'k--');
plot ([fs1 fs1]*fold*L,[-80 5],'k:');
% plot ([fg1 fg1],[-80 5],'k-');
axis ([0 fold*L/2 -80 5]);
xlabel('F [Hz]');
title('filtr 1');

for i=2:3,
    set(i, 'color', 'w');
end



N = nn+1;

end
function N = id_2(nr_sygnalu, fold, fnew, Li, Mi, FPi, FSi)

L = fnew / gcd(fold, fnew);
M = fold / gcd(fold, fnew);

L1 = Li(1);
M1 = Mi(1);
L2 = Li(2);
M2 = Mi(2);

fp1 = FPi(1);
fs1 = FSi(1);
fp2 = FPi(2);
fs2 = FSi(2);



% %definiowanie wejœciowej i wyjœciowej szybkoœci próbkowania
% fold = 320;
% fnew = 10240;
% 
% %obliczanie krotnoœci
% L = fnew/gcd(fold,fnew);
% M = fold/gcd(fold,fnew);
% L1 = 32;
% L2 = L/L1;
% M1 = 1;
% M2 = M/M1;


%definiowanie sygna³u testowego
x = gen_signal(nr_sygnalu, fold, L, M);

%obliczanie zafalowan w pasmach przepustowym i zaporowym
Ap = 0.1;
dp = -((1-10^(Ap/20))/(1+10^(Ap/20)));
As = 60;
ds = 10^(-As/20);

%uzupe³nianie zerami I stopien (interpolacja – wstawiamy L1-1 zer pomiedzy ka¿d¹ parê próbek)
xx = zeros(1,L1*length(x));
xx(1:L1:length(x)*L1)=x;

[nn1,fo1,mo1,w1] = remezord( [fp1 fs1], [1 0], [dp/2 ds], 1);
if nn1 < 3, nn1=3; end;
[nn,fo,mo,w] = remezord( [fp2 fs2], [1 0], [dp/2 ds], 1);
if nn < 3, nn=3; end;
N = [ nn1, nn ];
if sum(N) > 2000,
  for ind = 2:7
    subfig = findobj('tag', sprintf("MultiSRC_gui_v2b_fig%i", ind));
    if ~isempty(subfig)
      set(subfig, 'visible', 'off');
    end
  end

  errordlg('wymagany sumaryczny rz¹d filtrów > 2000', 'Zbyt du¿y rz¹d filtrów')
  return
end


%projektowanie filtru

h1 = L1*remez(nn1,fo1,mo1,w1);
% [H1, F21] = freqz(h1/L1, 1, 8192*32, 1);
[H1, F21] = freqz(h1/L1, 1, 8192*32, fold*L1);

%splot sygna³u po interpolacji z odpowiedzi¹ impulsow¹ filtru stopien I
y = conv(h1,xx);

% usuwanie stanu przejsciowego i decymacja sygba³u - wybieramy co M1-ta
% probke
p = length(h1)-1;
yy=y(p+1:M1:end);

%uzupe³nianie zerami stopien II
xx2 = zeros(1,L2*L1*length(x)/M1);
xx2(1:L2:length(x)*L2*L1/M1)=yy;

%projektowanie filtru stopien II
h2 = L2*remez(nn,fo,mo,w);
[H2, F2] = freqz(h2/L2, 1, 8192*32, fold*L1/M1*L2);

%splot sygna³u po interpolacji z odpowiedzi¹ impulsow¹ filtru stopien II
y2 = conv(h2,xx2);

% usuwanie stanu przejsciowego i decymacja sygna³u - wybieramy co M2-ta
% probke
p2 = length(h2)-1;
yy2=y2(p2+1:M2:end);

%wykresy
% figure(12)
% subplot(1,2,1);
% plot(x);
% subplot(1,2,2);
% plot(yy2);

%widmo sygna³u wejœciowego
wind0 = blackman(length(x));
Y0 = 20*log10(abs(fftshift(fft(wind0'.*x,32768))));
F_Y0 = linspace(-fold/2, fold/2, 32768+1); F_Y0(end) = [];
fig = findobj('Tag', 'MultiSRC_gui_v2b_fig2');
if ~isempty(fig)
  figure(fig);
  clf(fig);
else
  fig = figure('IntegerHandle', 'off', 'Tag', 'MultiSRC_gui_v2b_fig2', 'Color', 'w');
end
subplot(1,2,1);
plot(F_Y0, Y0);
set(gca, 'xlim', [-fold/2, fold/2]);
xlabel('F [Hz]');
title('syg. wejsciowy');

%widmo sygba³u wyjsciowego
figure(fig)
subplot(1,2,2);
wind1 = blackman(length(yy2));
Y1 = 20*log10(abs(fftshift(fft(wind1'.*yy2,32768))));
F_Y1 = linspace(-fnew/2, fnew/2, 32768+1); F_Y1(end) = [];
plot(F_Y1, Y1);
set(gca, 'xlim', [-fnew/2, fnew/2]);
xlabel('F [Hz]');
title('syg. wyjsciowy');

%filtr I
fig = findobj('Tag', 'MultiSRC_gui_v2b_fig3');
if ~isempty(fig)
  figure(fig);
  clf(fig);
else
  fig = figure('IntegerHandle', 'off', 'Tag', 'MultiSRC_gui_v2b_fig3', 'Color', 'w');
end
plot(F21,20*log10(abs(H1)));
hold on;
plot ([fp1 fp1]*fold*L1,[-80 5],'k--');
plot ([fs1 fs1]*fold*L1,[-80 5],'k:');
% plot ([fg1 fg1]*fold*L1,[-80 5],'k-');
axis ([0 0.5*fold*L1 -80 5]);
xlabel('F [Hz]');
title('filtr I');

%filtr II
fig = findobj('Tag', 'MultiSRC_gui_v2b_fig4');
if ~isempty(fig)
  figure(fig);
  clf(fig);
else
  fig = figure('IntegerHandle', 'off', 'Tag', 'MultiSRC_gui_v2b_fig2', 'Color', 'w');
end
plot(F2,20*log10(abs(H2)));
hold on;
plot ([fp2 fp2]*fold*L1/M1*L2,[-80 5],'k--');
plot ([fs2 fs2]*fold*L1/M1*L2,[-80 5],'k:');
% plot ([fg2 fg2]*fold*L1/M1*L2,[-80 5],'k-');
axis ([0 0.5*fold*L1/M1*L2 -80 5]);
xlabel('F [Hz]');
title('filtr II');

%Filtr zbiorczy
h2L = zeros(1,length(h2)*M1);
h2L(1:M1:length(h2)*M1)= h2;

h1L = zeros(1,length(h1)*L2);
h1L(1:L2:length(h1)*L2)= h1;

t=conv(h1L,h2L);

fig = findobj('Tag', 'MultiSRC_gui_v2b_fig5');
if ~isempty(fig)
  figure(fig);
  clf(fig);
else
  fig = figure('IntegerHandle', 'off', 'Tag', 'MultiSRC_gui_v2b_fig5', 'Color', 'w');
end
[H2L] = freqz(h2L/L2, 1, 8192*32, fold*L);
[H1L] = freqz(h1L/L1, 1, 8192*32, fold*L);
[H12, F212] = freqz(t/L, 1, 8192*32, fold*L);
subplot(2,1,1);
plot(F212,20*log10(abs(H1L)));
hold on
plot(F212,20*log10(abs(H2L)), 'r');
hold off
axis ([0 0.5*fold*L -80 5]);
xlabel('F [Hz]');
title('charakterystyki filtrów sk³adowych po przeniesieniu ich na szybkoœæ Fp1');
subplot(2,1,2);
plot(F212,20*log10(abs(H12)));
axis ([0 0.5*fold*L -80 5]);
xlabel('F [Hz]');
title('charakterystyka zbiorcza');


fig = findobj('Tag', 'MultiSRC_gui_v2b_fig6');
if ~isempty(fig)
  set(fig, 'visible', 'off');
end
fig = findobj('Tag', 'MultiSRC_gui_v2b_fig7');
if ~isempty(fig)
  set(fig, 'visible', 'off');
end

N = [ nn1+1, nn+1 ];

end
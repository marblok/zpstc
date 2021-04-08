function h=design_zad_2

%interpolation filter
Fp1 = 22050; % wej�ciowa szybko�� pr�bkowania
Fp2 = 44100; % wyj�ciowa szybko�� pr�bkowania
L = 2; % krotno�� interpolacji (Fp2/Fp1)

% zafalowanie w pasmie +/-0.1 dB, t�umienie poza pasmem -96dB
dp = 10.^(0.1/20)-1;
ds = 10.^(-96/20);

% Pasmo przej�ciowe 10-12kHz, filtr pracuj�cy na szybko�ci Fp2
c = firpmord( [10000 12000], [1 0], [dp ds], Fp2, 'cell');
h = firpm(c{:});
N= length(h)

figure(1)
plot(h)
% pause

figure(2)
freqz(h,1, 8*2048, L*Fp1)

% zapisanie wsp�czynnik�w filtru do pliku
dane.h  = L*h;
dane.Fp = Fp1;
save_filter_coef('cw1_zad2', dane);


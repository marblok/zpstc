function h=design_task_2

%interpolation filter
Fp1 = 22050; % input sampling rate
Fp2 = 44100; % output sampling rate
L = 2; % interpolation ratio (Fp2/Fp1)

% passband ripples +/-0.1 dB, out of band attenuation -96dB
dp = 10.^(0.1/20)-1;
ds = 10.^(-96/20);

% transition band 10-12kHz, filter works at rate equal to Fp2
c = firpmord( [10000 12000], [1 0], [dp ds], Fp2, 'cell');
h = firpm(c{:});
N= length(h)

figure(1)
plot(h)
% pause

figure(2)
freqz(h,1, 8*2048, L*Fp1)

% writing filter coefficients to file
dane.h  = L*h;
dane.Fp = Fp1;
save_filter_coef('ex1_task2', dane);


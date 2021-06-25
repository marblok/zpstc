function h=design_task_1

%interpolation filter
Fp1 = 22050; Fp2 = 44100;
L = 2;

% passband ripple +/-0.1 dB, out of band attenuation -96dB
dp = 10.^(0.1/20)-1;
ds = 10.^(-96/20);

c = firpmord( [10000 12000], [1 0], [dp ds], Fp2, 'cell');
h = firpm(c{:});
N= length(h)

figure(1)
plot(h)
% pause

figure(2)
freqz(h,1, 8*2048, L*Fp1)


dane.h  = L*h;
dane.Fp = Fp1;
save_filter_coef('ex2_task1', dane);


function h = design_task_2

% Design: interpolation filter
Fp_symb = 2000;
Fp_out = 8000;
L = Fp_out/Fp_symb

N = 10*L+1;
% filter type: square root raised cosine (podniesiony spierwiastkowany kosinus)
h = 0.5* 1/L*rcos4(N, L, 0.33);


figure(12)
set(12, 'unit', 'normalized', 'position', [0.05, 0, 1-0.1, 2/3-0.1])
freqz(h,1, 8*2048, L*Fp_symb)

figure(11)
set(11, 'unit', 'normalized', 'position', [0.05, 2/3, 1-0.1, 1/3-0.1])
plot(h)

dane.h  = L*h;
dane.Fp = Fp_symb;
save_filter_coef('ex2_task2', dane);


function ex5_task1

% Za³o¿enia na algorytm zmiany szybkoœci próbkowania
Fp1 = 44100; % wejœciowa szybkoœæ próbkowania
Fp2 = 48000; % wyjœciowa szybkoœæ próbkowania
F_max = 20000 % maksymalna czêstotliwoœæ pasma przenoszona wy³¹cznie ze zniekszta³ceniami liniowymi

AdB = 0.1; % Maksymalny poziom zniekszta³ceñ liniowych w pasmie do F_max
BdB = 96;  % T³umienie sk³adowych pochodz¹cych ze zniekszta³ceñ aliasowych w pasmie do F_max

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
M = Fp1 / gcd(Fp1,Fp2) 
L = Fp2 / gcd(Fp1,Fp2)

% Parametry filtru I-FIR 
L_ifir = 80; % wiekoœæ zeroinsertingu
Fp_sh = Fp1*L/L_ifir
Fp_ir = Fp1*L

F_c_prot = min(Fp1, Fp2)/2
F_c_ir = Fp_sh - min(Fp1, Fp2)/2

% projektowanie filtru prototypowego (kszta³tuj¹cego)
Ap = AdB/2; As = BdB;
[N_prot, h_prot] = get_filter(F_max, F_c_prot, Ap, As, Fp_sh, 1);
N_prot

% projektowanie filtru usuwaj¹cego repliki
[N_ir, h_ir] = get_filter(F_max, F_c_ir, Ap, As, Fp_ir, 1);
N_ir


figure(1);
set(1, 'unit', 'normalized' , 'position', [0.05, 0.05, 1-0.1, 1-0.15]);
subplot(2,2,1);
plot(h_prot);
title('odpowiedŸ impulsowa filtru prototypowego');
subplot(2,2,2);
[H, F] = freqz(h_prot,1, 8*2048, L/L_ifir*Fp1);
plot(F, 20*log10(abs(H)));
title('charakterystyka czêstotliwoœciowa filtru prototypowego');

subplot(2,2,3);
plot(h_ir);
title('odpowiedŸ impulsowa filtru usuwaj¹cego repliki');
subplot(2,2,4);
[H, F] = freqz(h_ir,1, 8*2048, L*Fp1);
plot(F, 20*log10(abs(H)));
title('charakterystyka czêstotliwoœciowa filtru usuwaj¹cego repliki');

dane.h{1}  = h_prot;
dane.h{2}  = L_ifir;
dane.Fp = L/L_ifir*Fp1;
save_filter_coef('ex5_task1_h_sh', dane, 1);

% format long
% h_sh(1:13)
% return
dane.h  = h_ir;
dane.Fp = L*Fp1;
save_filter_coef('ex5_task1_h_ir', dane, 1);

% wyznaczenie odpowiedzi impulsowej filtru zbiorczego
h_sh = zeros(1,length(h_prot)*L_ifir);
h_sh(1:L_ifir:end) = h_prot;

h_all = conv(h_sh, h_ir);
N_all = length(h_all)

figure(2);
set(2, 'unit', 'normalized' , 'position', [0.05, 0.05, 1-0.1, 1-0.15]);
subplot(2,2,1);
plot(h_all);
title('odpowiedŸ impulsowa filtru zbiorczego');
subplot(2,2,2);
[H, F] = freqz(h_all,1, 8*2048, L*Fp1);
plot(F, 20*log10(abs(H)));
title('charakterystyka czêstotliwoœciowa filtru zbiorczego');

% h_all_2 = fileread('ex5_task1.flt');
h_all_2 = fileread('ex5_task1.wav');
subplot(2,2,3);
plot(h_all_2);
title('odpowiedŸ impulsowa filtru zbiorczego');
subplot(2,2,4);
[H, F] = freqz(h_all_2,1, 8*2048, L*Fp1);
plot(F, 20*log10(abs(H)));
title('charakterystyka czêstotliwoœciowa filtru zbiorczego');

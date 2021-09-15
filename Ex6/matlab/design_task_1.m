function h=design_task_1
% type = 2;
% 
% %interpolation filter
% switch type,
%   case 1,
%     Fp1 = 22050; Fp2 = 44100;
%     
%   case 2,
%     Fp1 = 8000; Fp2 = 32000;
%     F_trans = [3000 5000];
% end
% L = Fp2 / Fp1;

% M = 1;
Fsymb = 1500; Fp_out = 32*Fsymb;

L = Fp_out/Fsymb


% filtr dla pierwszego stopnia interpolacji
r = 0.33;
h_rc=rcos4(10*L,L, 0.33);

K = 8192;
[H_rc, f] = freqz(h_rc/L, 1, K, L*Fsymb);

figure(1)
subplot(3,1,1)
plot(h_rc)
subplot(3,1,2)
plot(f, abs(H_rc))
set(gca, 'Ylim', [0, 1.1], 'Xlim', [0, L*Fsymb/2])
subplot(3,1,3)
plot(f, 20*log10(abs(H_rc)))
set(gca, 'Ylim', [-70, 3], 'Xlim', [0, L*Fsymb/2])

dane.h  = h_rc/L;
dane.Fp = Fsymb;
save_filter_coef('ex6_task1_h_rc', dane, 1);



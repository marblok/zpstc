function h=design_task_3
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

Fsymb = 2400; Fp_out = 48000;

L = Fp_out/Fsymb

L1 = 4
L2 = L / L1

% filtr dla pierwszego stopnia interpolacji
r = 0.33;
h_rc=rcos4(20*L1,L1, 0.33);

K = 8192;
[H_rc, f] = freqz(h_rc/L1, 1, K, L1*Fsymb);

figure(1)
subplot(3,1,1)
plot(h_rc)
subplot(3,1,2)
plot(f, abs(H_rc))
set(gca, 'Ylim', [0, 1.1], 'Xlim', [0, L1*Fsymb/2])
subplot(3,1,3)
plot(f, 20*log10(abs(H_rc)))
set(gca, 'Ylim', [-70, 3], 'Xlim', [0, L1*Fsymb/2])

% filtr dla drugiego stopnia interpolacji
Fd = Fsymb/2 * (1+0.33);
% !!! pasmo przepustowe filtru drugiego stopnia obejmuj�ce 
% pasmo przepustowe i przej�ciowe filtru pierwszego stopnia
% !!! jest to nietypowe rozwi�zanie gwarantuj�ce zachowanie charakteru
% pasma przej�ciowego filtru kszta�tuj�cego
Fg = L1*Fsymb - Fd;

% zafalowanie w pasmie +/-0.1 dB, t�umienie poza pasmem -60dB
dp = 10.^(0.1/20)-1;
ds = 10.^(-60/20);

c = firpmord( [Fd, Fg], [1 0], [dp ds], Fp_out, 'cell');
h2 = firpm(c{:});
N2 = length(h2)
[H2, f] = freqz(h2, 1, K, Fp_out);

figure(2)
subplot(3,1,1)
plot(h2)
subplot(3,1,2)
plot(f, abs(H2))
set(gca, 'Ylim', [0, 1.1], 'Xlim', [0, Fp_out/2])
subplot(3,1,3)
plot(f, 20*log10(abs(H2)))
set(gca, 'Ylim', [-70, 3], 'Xlim', [0, Fp_out/2])

% analiza filtru zbiorczego
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% filtr h_rc - pracuje na L1*Fsymb
% przenosimy go na szybko�� pr�bkowania filtru drugiego
h1 = zeros(1, length(h_rc)*L2);
h1(1:L2:end) = h_rc;
[H1, f] = freqz(h1/L1, 1, K, Fp_out);

figure(3)
subplot(2,1,1)
plot(f, abs(H1), 'r');
hold on
plot(f, abs(H2), 'b');
plot(f, abs(H1.*H2), 'k');
hold off
set(gca, 'Ylim', [0, 1.1], 'Xlim', [0, Fp_out/2])
subplot(2,1,2)
plot(f, 20*log10(abs(H1)), 'r');
hold on
plot(f, 20*log10(abs(H2)), 'b');
plot(f, 20*log10(abs(H1.*H2)), 'k'); % filtr zbiorczy
hold off
set(gca, 'Ylim', [-80, 3], 'Xlim', [0, Fp_out/2])

dane.h  = h_rc;
dane.Fp = L1*Fsymb;
save_filter_coef('cw3_zad3_h_rc', dane, 1);


dane.h  = L2*h2;
dane.Fp = Fp_out;
save_filter_coef('cw3_zad3_h2', dane, 1);


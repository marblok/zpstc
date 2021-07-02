function gen_sign_task1

% delta_44100.wav
Fp = 44100; 
x = zeros(1*Fp,1); x(1) = 1;
try 
  wavwrite(x, Fp, 16, 'delta_44100.wav');
catch
  audiowrite('delta_44100.wav', x, Fp, 'BitsPerSample', 16);
end



% % test1_44100.wav, test2_44100.wav
% Fp = 44100; Fu = 0.4*Fp;
% x = gen_signal(Fu, 1, Fp);
% wavwrite(x, Fp, 16, 'test1_44100.wav');
% 
% Fp = 44100; Fu = 0.5*Fp;
% x = gen_signal(Fu, 1, Fp);
% wavwrite(x, Fp, 16, 'test2_44100.wav');

function x = gen_signal(Fu, T, Fp)

N = T*Fp;
n = 0:N-1;
f = linspace(0, Fu, N);
faza = 2*pi*cumsum(f)/Fp;

x = 0.9*cos(faza);


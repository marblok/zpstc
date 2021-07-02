function gen_sign_task2

% % test_signal.wav
% Fp = 8000; Fu = 0.5*Fp;
% x = gen_signal(Fu, 1, Fp);
% wavwrite(x, Fp, 16, 'test_signal.wav');
% 
% % test1_11025.wav, test2_11025.wav
% Fp = 11025; Fu = 0.4*Fp;
% x = gen_signal(Fu, 1, Fp);
% wavwrite(x, Fp, 16, 'test1_11025.wav');
% 
% Fp = 11025; Fu = 0.5*Fp;
% x = gen_signal(Fu, 1, Fp);
% wavwrite(x, Fp, 16, 'test2_11025.wav');


% test_44100.wav, test_48000.wav
Fp = 44100; Fu = 0.5*Fp;
x = gen_signal(Fu, 1, Fp);
try 
  wavwrite(x, Fp, 16, 'test_44100.wav');
catch
  audiowrite('test_44100.wav', x, Fp, 'BitsPerSample', 16);
end

Fp = 48000; Fu = 0.5*Fp;
x = gen_signal(Fu, 1, Fp);
try
  wavwrite(x, Fp, 16, 'test_48000.wav');
catch 
  audiowrite('test_48000.wav', x, Fp, 'BitsPerSample', 16);
end

function x = gen_signal(Fu, T, Fp)

N = T*Fp;
n = 0:N-1;
f = linspace(0, Fu, N);
faza = 2*pi*cumsum(f)/Fp;

x = 0.9*cos(faza);


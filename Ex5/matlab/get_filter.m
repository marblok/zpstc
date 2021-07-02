function [N, h] = get_filter(F_p, F_s, Ap, As, Fp, tryb)
% tryb opcjonalnie
%   tryb = 0 tylko szacowanie d³ugoœci odpowiedzi impulsowej
%   tryb = 1 szacowanie d³ugoœci odpowiedzi impulsowej i projektowanie filtru
if nargin == 5,
  tryb = 0;
end

% zafalowanie w pasmie +/-Ap dB, t³umienie w paœmie zaporowm 96dB
dp = 10.^(Ap/20)-1; ds = 10.^(-As/20);

% szacowanie parametrów filtru
c = firpmord( [F_p, F_s], [1 0], [dp ds], Fp, 'cell');

if tryb == 0,
  h = [];
  N = c{1};
  return;
end

K = 8*2048;

dN = 1; old_PE = 0;
while(1),
  h = firpm(c{:});
  [H, F] = freqz(h,1, K, Fp);
  H_dB = 20*log10(abs(H(1:K/2+1)));
  ind = find(F(1:K/2+1) >= F_s);
  PE = max(H_dB(ind));
  
  if (PE <= -As)
    break
  end
  if old_PE < 0,
    dPE = PE - old_PE;
    dN = round((-As - PE)/dPE);
  end
  old_PE = PE;
  c{1} = c{1} + dN;
end

N = length(h);

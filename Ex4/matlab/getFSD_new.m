% %poszukiwanie filtru optymalnego o dowolnej d³ugoœci
%w oparciu o uk³ad równañ rzeczywistych
% w oparciu o phase rotated real error
%
%implementacja w oparciu o poszukiwanie ekstremów
%pomiêdzy dwoma kolejnymi punktami ekstremalnymi o tym samym znaku
% 
% Wariant bezrastrowy w u¿yciem interpolacji parabolicznej
%
%getFSD_new(N,voff,d); 
% N - d³ugoœæ odpowiedzi impulsowej
% voff - czêstotliwoœæ unormowana koñca pasma dobrej aproksymacji
% d - opóŸnienie netto

% 2011.03.01
function [h, TPEopt, ekst, A, B]=getFSD_new(N, voff, d, K, ekst, force_end);

if nargin < 6,
  force_end = 0;
end

if nargin==0,
  N=4;
  voff=0.35;
  d=0.25;
  K=2048; %DFT length
elseif nargin==1,
  voff=0.35;
  d=0.25;
  K=2048; %DFT length
elseif nargin==2
  d=0.25;
  K=2048; %DFT length
elseif nargin==3
  K=2048; %DFT length
end;

% df_factor = 10; df_factor_scaling = 1.05;
df_factor = 20; df_factor_scaling = 1.0;
df_factor = 4;

% % K = K * 8 * 16;
% K_ = 8192*2*64;

%   N_LPF(ind)=ceil((-20*log10(sqrt(dp*ds))-13)/(14.6*(vs-vp))+1);
% N_FSD_est=ceil(-1/8*(PE/(4*(1/2-voff))+3));
PE_est=-(N*8+3)*(4*(1/2-voff));
if (PE_est < -140)
  warning('Peek error estimation based on impulse response length shows that N is so large that it might result in fatal design errors');
elseif (PE_est < -200)
  warning('Peek error estimation based on impulse response length shows that N is so large that it probably will result in fatal design errors');
elseif (PE_est < -275)
  warning('Peek error estimation based on impulse response length shows that N is so large that it must result in fatal design errors');
  pause
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bez kwantyzacji osi czêstotliwoœci
M_in=round(voff*K)+1; % voff=(M_in-1)/K;

L=(N-1)/2; % bulk delay
H_L = Hideal(L, voff, K); H_L=conj(H_L);
Hid = Hideal(L+d, voff, K);

TPEopt=-300;
if nargin<5
% 	if rem(N,2)==0
% 	  v=linspace(0,voff,N/2+1);
% 	else
% 	  v=linspace(voff/N,voff,(N+1)/2);
% 	end;
%  alfa=0.95;
%  alfa=0.975;
  alfa=0.999;
  if (voff<=0.47)
    alfa=0.99;
  end
%  alfa=0.996;
%  alfa=0.999;
	if rem(N,2)==0
%       '0'
    x=[0:(N/2-1)];
    v=[0 cumsum(alfa.^x)];
    v=v/v(end)*voff;
	else
    x=[0:((N-1)/2)];
    v=cumsum(alfa.^x);
    v(1)=v(1)/2;
    v=v/v(end)*voff;
	end;
  ekst=(v*K);
%   ekst=round(v*K);
end

nn=1:N; A=zeros(N+1,N+1); B=zeros(N+1,1);
if (rem(N,2)==0) %parzyste
  N_2=N/2;
  A(1:N_2,N+1)=zeros(N_2,1);
  A(((N_2+1):N),N+1)=((-1).^(1:N_2).');
  A(N+1,(1:(N+1)))=ones(1,N+1);
  B(N+1)=1;
else
  N_2=(N+1)/2;
  A(1:N_2,N+1)=-((-1).^(1:N_2).');
  A(((N_2+1):(N+1)),N+1)=0;
end

% ile = 0;
koniec = 1;
while 1,
%   ile = ile + 1
  %construct equations  
  if (rem(N,2)==0) %parzyste
    v=ekst(2:length(ekst))/K; %N-parzyste
    for k=1:N_2,
      A(k,nn)=sin(2*pi*(nn-1-L)*v(k));
      B(k)=sin(2*pi*v(k)*d);
      A((N_2+k),nn)=cos(2*pi*(nn-1-L)*v(k));
      B((N_2+k))=cos(2*pi*v(k)*d);
    end
  else
    v=ekst(1:length(ekst))/K; %N-nieparzyste
    for k=1:N_2,
      A(k,nn)=sin(2*pi*(nn-1-L)*v(k));
      B(k)=sin(2*pi*v(k)*d);
      A((N_2+k),nn)=cos(2*pi*(nn-1-L)*v(k));
      B((N_2+k))=cos(2*pi*v(k)*d);
    end
  end
  
  if any(~isfinite(B)),
    disp('B isn''t finite');
    error
  end
  if any(~isfinite(A)),
    disp('A isn''t finite');
    error
  end
  x=A\B; %  x=inv(A)*B;
  h=x(1:N).';
  if any(~isfinite(h)),
    disp('h isn''t finite');
    plot(20*log10(abs(Error)));
    error
  end
  
  if any(~isfinite(h)),
    disp('Impulse response isn''t finite');
    error
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Error calculation
  H = fft(h,K);
  H = H(1:M_in);
  H_L = H_L(1:M_in);
  Hid = Hid(1:M_in);

  if (M_in-1)/K < voff,
    n = 0:length(h)-1;
    H_off = sum(h.*exp(-i*2*pi*voff*n));
    Hid_off = exp(-i*2*pi*voff*(L+d));
    H_L_off = exp(+i*2*pi*voff*L);
    H = [H, H_off];
    Hid = [Hid, Hid_off];
    H_L = [H_L, H_L_off];
  end
  
  Error=H-Hid;
  Error=Error.*H_L;

  n = 0:length(h)-1;
  H_off = sum(h.*exp(i*2*pi*voff*(n-L).*n));
  Hid_off = exp(-i*2*pi*voff*d);
  
  % wyznaczenie wspó³czynnika rotacji b³êdu
  % => korygujemy b³¹d tak, ¿eby na czêstotliwoœci 
  Error_v_off = H_off - Hid_off;
%   if abs(Error(M_in)) > 0,
  if abs(Error_v_off) > 0,
    Err_rot = conj(Error_v_off)/abs(Error_v_off);
%   Err_rot = conj(Error(M_in))/abs(Error(M_in));
    Error=Error.*Err_rot;
  else
    Err_rot = 1;
  end
  if any(~isfinite(Error)),
    disp('Normalized error isn''t finite');
    error
  end

  if force_end == 1;
    TPEopt = 20*log10(abs(x(N+1)));
    TPEopt(2) = 20*log10(max(abs(Error)));
    return;
  end
%   plot(20*log10(abs(Error)));
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % poszukiwanie po³o¿enia nowych ekstremów
%   ekst=Ekstr(real(Error),ekst, N);
%  ekst=Ekstr_new(h, ekst, Err_rot, d, N);
  ekst_old = ekst;
  ekst=Ekstr_new(h, ekst/K, d, N, K, df_factor)*K;
  df_factor = round(df_factor * df_factor_scaling);
  if df_factor > 200,
    df_factor = 200;
  end
  
  if any(~isfinite(ekst)),
    disp('ekst isn''t finite');
    error
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % wyznaczenie peek error
  %TPE=max(20*log10(abs(Error(ekst+1)))); 
  %TPE_H=max(20*log10(abs(Error))); 
  TPE=20*log10(abs(x(N+1)));
  if TPEopt<TPE
    TPEopt=TPE;
  else
%     ekst
    ekst_old - ekst;
    koniec = koniec - 1;
    if koniec == 0,
      break;
    end
  end
end

return


%========================================================
function Hid = Hideal(delay, F, K)
% delay - opóŸnienie u³amkowe projektowanego filtru
% F - unormowana czêstotliwoœæ górna pasma akuratnej aproksymacji
% K - liczba pr¹¿ków transformaty
% Hid - charakterystyka idealna (K-punktowa)

v=(0:(K-1))/K; v(K/2:K)=v(K/2:K)-1;
Hid = ones(1,K).*exp(-i*2*pi*v*delay);

%========================================================
function x=Ekstr(Err, old, N)
% je¿eli filtr parzystej d³ugoœci to zak³adamy, ¿e punkty ekstremalne w 0 i M_in (end)
% nie ulegaj¹ zmianie
% je¿eli filtr o nieparzystej d³ugoœci to w zerze napewno nie ma punktu ekstremalnego 
% oraz punkt ekstremalny w M_in (end) nie ulega zmianie

delta=max(abs(Err(old+1)));

%wyjœciowy zestaw punktów ekstremalnych
temp=old+1; %parzyste
if rem(N,2)==1
  temp=[1 temp]; %nieparzyste
end

min_max=0;
for ind=length(temp):-1:3
  pom=Err(temp(ind-2):temp(ind));
  if min_max==0
    [wart, ind_new]=find(pom==min(pom));
  else
    [wart, ind_new]=find(pom==max(pom));
  end
  if abs(wart)>delta
    temp(ind-1)=temp(ind-2)+ind_new-1;
  end
  min_max=1-min_max;
end

x=temp-1;
if rem(N,2)==1
	x=x(2:length(x));
end

function f_new = Ekstr_new(h, f_old, d, N, K, df_factor);
% je¿eli filtr parzystej d³ugoœci to zak³adamy, ¿e punkty ekstremalne w 0 i M_in (end)
% nie ulegaj¹ zmianie
% je¿eli filtr o nieparzystej d³ugoœci to w zerze napewno nie ma punktu ekstremalnego 
% oraz punkt ekstremalny w M_in (end) nie ulega zmianie
L = (N-1)/2;
n = 0:N-1; n=n(:);

f_old_ = f_old;
foff = f_old(end);
H_off = sum(h(:).*exp(-i*2*pi*foff*(n-L)));
Hid_off = exp(-i*2*pi*foff*d);
Error_f_off = H_off - Hid_off;
if abs(Error_f_off) > 0,
  Err_rot = conj(Error_f_off)/abs(Error_f_off);
else
  Err_rot = 1;
end

for ind = 1:length(f_old),
  H_off(ind) = sum(h(:).*exp(j*2*pi*f_old(ind)*(n-L).*n));
end
Hid_old = exp(-i*2*pi*f_old*d);

Err_old = real((H_off - Hid_old)*Err_rot);

% delta=max(abs(Err(old+1)));
delta=max(abs(Err_old));


do_draw = 0;
if do_draw == 1,
  % %   f = linspace(0,0.5, 2048);
  f = linspace(0,0.5, K); % 2011.02.28
  
  E_ = 0;
  H_ = 0;
  for n_ = n.',
    E_ = E_ + Err_rot*(-j)*2*pi* h(n_+1).*(n_-(L+d)).*exp(-j*2*pi*f*(n_-(L+d)));
    H_ = H_ + Err_rot*h(n_+1).*exp(-j*2*pi*f*(n_-L));
  end

  H_id = Err_rot*exp(-j*2*pi*f*d);

  figure(5)
  subplot(3,1,1)
  plot(f, real(E_));
  hold on
  plot(f, imag(E_), 'r');
  hold off
  set(gca, 'Xtick', f_old_, 'Xgrid', 'on');
  subplot(3,1,2)
  plot(f, abs(H_-H_id), 'k');
  hold on
  plot(f, real(H_-H_id), 'b');
  plot(f, imag(H_-H_id), 'r');
  hold off

  subplot(3,1,3)
  plot(f, 20*log10(abs(H_-H_id)), 'k');
end
  
if rem(N, 2) == 1,
  f_old = [0, f_old]; % 0 must always be included 
end 
% find extremas
% leave point at v_off (maximum)
f_new = f_old(end);

polarization = -1;
while length(f_old) > 2,
  % find extremum between f_start and f_end 
  f_end = f_old(end);
  f_start = f_old(end-2);

  fo = f_old(end-1); % old value
  
  df = min([(f_end-fo), (fo-f_start)]) / df_factor;
  
  fo_new = find_extremum(h, d, Err_rot, [fo-df, fo, fo+df], polarization, do_draw);
  
  f_new = [fo_new, f_new];
  f_old(end) = [];
  
  polarization = -polarization;
end
if rem(N, 2) == 0,
  f_new = [f_old(1), f_new]; % include 0
end

% size(f_old_)
% size(f_new)

if do_draw == 1,
  subplot(3,1,2)
  set(gca, 'Xtick', f_new, 'Xgrid', 'on');
  subplot(3,1,3)
  set(gca, 'Xtick', f_new, 'Xgrid', 'on');

  % f_new
  % 
  pause
end

function fo_new = find_extremum(h, d, Err_rot, fs, polarization, do_draw)

N = length(h);
n = 0:N-1; n = n(:);
L = (N-1)/2;

f = fs;
% f = [fs(1), (fs(1)+fs(2))/2, fs(2), (fs(2)+fs(3))/2, fs(3)];
% f = linspace(fs(1), fs(3), 11);

% E_ = 0;
H_ = 0;
for n_ = n.',
%   E_ = E_ + Err_rot*(-j)*2*pi* h(n_+1).*(n_-(L+d)).*exp(-j*2*pi*f*(n_-(L+d)));
  H_ = H_ + Err_rot*h(n_+1).*exp(-j*2*pi*f*(n_-L));
end
H_id = Err_rot*exp(-j*2*pi*f*d);

E_ = real(H_-H_id) * polarization;

%f(x) = a*x^2+b*x+c
%delta = b^2 - 4*a*c
%f_ekstremum = -b/(2*a)

% % y = polyval(a, x)
% y1 = a*x1^2 + b*x1 + c
% y2 = a*x2^2 + b*x2 + c
% y3 = a*x3^2 + b*x3 + c

a = polyfit(f, E_, 2);

fo_new = -a(2)/(2*a(1));

% 
% if (fo_new < fs(1)) 
%   fo_new = fs(1);
% end
% if (fo_new > fs(3)) 
%   fo_new = fs(3);
% end

if (fo_new < fs(1)) || (fo_new > fs(3))
  ind = find(max(E_) == E_);
  fo_new = f(ind(1));
else
  H_new = 0;
  for n_ = n.',
    H_new = H_new + Err_rot*h(n_+1).*exp(-j*2*pi*fo_new*(n_-L));
  end
  H_id_new = Err_rot*exp(-j*2*pi*fo_new*d);
  E_new = real(H_new-H_id_new) * polarization;

  E_ = [E_, E_new];
  f = [f, fo_new];
  
  ind = find(max(E_) == E_);
  if (fo_new ~= f(ind(1))),
    fo_new = f(ind(1));
  end

end

% % if fo_new < fs(1)
% %   fo_new = fs(1);
% % end
% % if fo_new > fs(3)
% %   fo_new = fs(3);
% % end
% % % fo_new
% % % fo_new;

% do_draw = 1;
if do_draw == 1,
  f_ = linspace(fs(1), fs(3), 51);
  y = polarization * polyval(a,  f_);
  subplot(3,1,2)
  hold on
  plot(f_, y, 'm');
  hold off
  max_y = max(abs(y));
  if max_y == 0, max_y = 1, end;
  set(gca, 'Ylim', [-1.15, 1.15]*max_y);
  % % pause
  % 
end

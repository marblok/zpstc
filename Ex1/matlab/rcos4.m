function h_rc=rcos4(N,L, r)
%square root raised cosine

fp=1000;

if nargin==2,
  r=0.5;
end

n=(0:N-1)-(N-1)/2;
% t=n/fp; T=L/fp;

% % L=pi*L; h_rc=sinc(pi*n/L).*cos(pi*r*n/L)./(L-4*(r.^2)*(n.^2)/L);
% h_rc=sinc(t/T).*cos(r*pi*t/T)./(1-4*(r.^2)*(t.^2)/T.^2);
h_rc=(1-r)*sinc((1-r)*n/L) + ...
     r*( sinc(r*n/L+1/4).*cos(pi*(n/L+1/4)) + ...
         sinc(r*n/L-1/4).*cos(pi*(n/L-1/4)) );
% % h_rc= 4*r* ( cos((1+r)*pi*n/L) + sin((1-r)*pi*n/L)./(4*r*n/L)) ...
% %       ./ (pi*sqrt(1/L)*(1-(4*r*n/L).^2));
% 
ind=~isfinite(h_rc);
h_rc(ind)=0;

if nargout==0,
	figure(1)
	subplot(3,1,1)
	stem(n/L,h_rc)
	subplot(3,1,2)
  n=-(N-1):(N-1);
%   size(n)
%   size(conv(h_rc,h_rc))
	stem(n/L,conv(h_rc,h_rc))
	subplot(3,1,3)
	H=freqz(h_rc,1,2048, 'whole');
	H=abs(H(1:end/2));
	f=linspace(0,1000,length(H));
	plot(f,H)
  
	H=freqz(ones(1,L),1,2048, 'whole');
	H=abs(H(1:end/2));
	f=linspace(0,1000,length(H));
  hold on
	plot(f,H,'k')
  hold off
	H=freqz(conv(h_rc,h_rc)/L,1,2048, 'whole');
	H=abs(H(1:end/2));
	f=linspace(0,1000,length(H));
  hold on
	plot(f,H,'m')
  hold off
end
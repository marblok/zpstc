function [h, PE_FSD] = FSDfilter(method, N, delay, voff_typ)
% delay - total filter delay

switch method
  case 'LSE'
    h=LSE(N,delay);
  case 'LSE2'
    h=LSE(N,delay,voff_typ);
  case 'Lagr'
    h=Lagr(N,delay);
  case 'window'
    h=window(N,delay,voff_typ);
  case 'chebwin'
    h=window(N,delay,3,voff_typ);
  case 'offset'
    h=offset(N,delay,voff_typ);
  case 'optimal'
    eval('h=getFSD(N,voff_typ,delay-(N-1)/2);','h=NaN;');
    
    if isnan(h)
      Krok=0;
      Factor=4;
      
      while isnan(h)
        Krok=Krok+1;
        krok=Krok; voff=voff_typ;
        while krok
          voff=voff+(1/2-voff)/Factor;
          krok=krok-1;
        end;
        eval('[h, PE_FSD, ekst]=getFSD(N,voff,delay-(N-1)/2);','h=NaN;');
        h=h.';
      end;
      
      while Krok>=0
        krok=Krok; voff=voff_typ;
        while krok
          voff=voff+(1/2-voff)/Factor;
          krok=krok-1;
        end;
        ekst=ekst/ekst(end)*voff;
        [h, PE_FSD, ekst]=getFSD(N,voff,delay-(N-1)/2, ekst);
        h=h.';
        Krok=Krok-1;
      end;
    end;
    
  case 'optimal2'
    % rasterless implementation
    [h, PE_FSD] =getFSD_new(N,voff_typ,delay-(N-1)/2);
%     eval('h=getFSD_new(N,voff_typ,delay-(N-1)/2);','h=NaN;');
    
    if isnan(h)
      Krok=0;
      Factor=4;
      
      while isnan(h)
        Krok=Krok+1;
        krok=Krok; voff=voff_typ;
        while krok
          voff=voff+(1/2-voff)/Factor;
          krok=krok-1;
        end;
        eval('[h, PE_FSD, ekst]=getFSD_new(N,voff,delay-(N-1)/2);','h=NaN;');
        h=h.';
      end;
      
      while Krok>=0
        krok=Krok; voff=voff_typ;
        while krok
          voff=voff+(1/2-voff)/Factor;
          krok=krok-1;
        end;
        ekst=ekst/ekst(end)*voff;
        [h, PE_FSD, ekst]=getFSD_new(N,voff,delay-(N-1)/2, ekst);
        h=h.';
        Krok=Krok-1;
      end;
    end;
end;

function h=Lagr(N,delay)

for n=0:N-1
  m=0:N-1;
  m(n+1)=[];
  h(n+1)=prod((delay-m)./(n-m));
end;

function h=LSE(N,delay, voff)

n=0:N-1;
if nargin ==2,
  h=sinc(n-delay);
elseif nargin==3
  alfa=2*voff;
  for k=0:N-1
    l=0:N-1;
    P(k+1,l+1)=alfa*sinc(alfa*(k-l));
    p(k+1,1)=alfa*sinc(alfa*(k-delay));
  end;
  h=inv(P)*p;
  h = h(:).';
end;

function h=window(N,delay,typ, Rp)
n=0:N-1;
hsinc=sinc(n-delay);

if typ==3
  wsym=chebwin(N,Rp).';
else
  switch typ
    case 2
      C=[0.43 -0.5  0.08]; %Blackman
    case 1
      C=[0.54 -0.46 0]; %Hamming
    case 0
      C=[0.5  -0.5  0]; %Hann
  end;
  wsym=C(1)+C(2)*cos(2*pi*n/(N-1))+C(3)*cos(4*pi*n/(N-1));
end;

h=hsinc.*wsym;

function h=offset(N,delay,typ)
n=0:N-1;
hsinc=sinc(n-delay);
switch typ
  case 2
    C=[0.43 -0.5  0.08]; %Blackman
  case 1
    C=[0.54 -0.46 0]; %Hamming
  case 0
    C=[0.5  -0.5  0]; %Hann
end;

epsilon=delay-(N-1)/2;
%wsym=C(1)+C(2)*cos(2*pi*n/(N-1))+C(3)*cos(4*pi*n/(N-1));
woffset=C(1)+C(2)*cos(2*pi*(n-epsilon)/(N-1))+C(3)*cos(4*pi*(n-epsilon)/(N-1));
h=hsinc.*woffset;


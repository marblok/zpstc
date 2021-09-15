function hf = eyediagram2(x, offset, Nsymb, hf)

if nargin == 4,
  figure(hf);
else
  hf = figure;
end
subplot(1,1,1)
if nargin == 2,
  Nsymb = offset;
  offset = ceil((Nsymb-1)/2);
%   Nsymb = length(x)+1;
end
x(ceil(length(x)/Nsymb)*Nsymb) = 0;
hp = [];
% t_ = linspace(-Nsymb/2, Nsymb/2, Nsymb);
t_ = 0:Nsymb;
ind_ = 1:(Nsymb+1);
% if rem(Nsymb,2) == 1,
%   t_ = (0:Nsymb+1) - Nsymb/2  - 0.5;
%   
%   ind_ = [0:Nsymb+1];
% else
%   t_ = (-1:Nsymb) - Nsymb/2 + 0.5;
%   
%   ind_ = [-1:Nsymb];
% end


dt = offset - round(offset);
t_ = t_ + dt; offset = round(offset);

% t_ = t_ / Nsymb;

if sum(abs(imag(x))) == 0,
  subplot(1,1,1);
  plot(t_, t_*0, 'k:');
  hold on
  plot([0, 0], [-1.1, 1.1], 'k:');
  % plot(t_, +1+t_*0, 'k:');
  % plot(t_, -1+t_*0, 'k:');
  for ind = offset:Nsymb:length(x)-(Nsymb+1),
    x_ = x(ind+ind_);
  
    hp(end+1) = plot(t_, real(x_));
%     hp(end+1) = plot(t_, imag(x_), 'r');
  end
  hold off
else
  subplot(1,2,1);
  plot(t_, t_*0, 'k:');
  hold on
  plot([0, 0], [-1.1, 1.1], 'k:');
  % plot(t_, +1+t_*0, 'k:');
  % plot(t_, -1+t_*0, 'k:');
  subplot(1,2,2);
  plot(t_, t_*0, 'k:');
  hold on
  plot([0, 0], [-1.1, 1.1], 'k:');
  for ind = offset:Nsymb:length(x)-(Nsymb+1),
    x_ = x(ind+ind_);
  
    subplot(1,2,1);
    hp(end+1) = plot(t_, real(x_));
    subplot(1,2,2);
    hp(end+1) = plot(t_, imag(x_), 'r');
  end
  hold off
end

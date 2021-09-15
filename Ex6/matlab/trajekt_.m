function klawisz=trajekt_(x_s, time, punkty, no_head, pause_on_start)
% \added 2019.10.18 ¿ó³te t³o okna w trybie pauzy
% \added 2019.10.18 prosty tryb testowania (ustawiaj¹c test_mode = 1; w linii 11)
% \fixed 2019.10.18 wy³¹czenie zoom przed definiowaniem 'KeyPressFcn'
% \added 2017.10.24 dodano obs³ugê trybu, w którym x_s zawiera dane i
%    etykiety (interpreter LaTeX): x_s{1}.data = ..., x_s{1}.label = '?';
% \fixed 2015.10.27 addapted to matlab 2015b
% \fixed 2006.11.14 fixed axes position in single and double plot mode
% \added 2006.11.14 input option no_head
% \added 2006.11.14 input option pause_on_start

test_mode = 0;
if (test_mode == 1) & (nargin == 0)
  x_s = exp(sqrt(-1)*0.1*[0:1000]);
  klawisz = trajekt_(x_s, 100, [], 0, 1);
  return;
end

if nargin < 5,
  pause_on_start = 0;
end
if nargin < 4,
  no_head = 1;
end
if nargin < 3,
  punkty = [];
end
labels = {};


if nargin == 0, 
error('Not enough input arguments.'); 
else
	if isreal(x_s),
    ile = 1;
    x = x_s;
    y = x; x = 1:length(y);
    
    clear x_s;
    x_s{1}.x = x;
    x_s{1}.y = y;
    max1=max(abs(x(~isinf(x))));
    zakres=[-max1 max1 -max1 max1];
	else
    if iscell(x_s) == 0,
      ile = 1;
      y=imag(x_s);
      x=real(x_s);
      clear x_s;
      x_s{1}.x = x;
      x_s{1}.y = y;
      max1=max(max(abs(x(~isinf(x)))),max(abs(y(~isinf(y)))));
      zakres=[-max1 max1 -max1 max1];
    else
      x_s_tmp = x_s;
      clear x_s
      
      ile = length(x_s_tmp);
      max1 = eps;
      for ind = 1:ile,
        if isfield(x_s_tmp{ind}, 'label'),
          labels{ind} = x_s_tmp{ind}.label;
        end
        if isfield(x_s_tmp{ind}, 'data'),
          x_s_tmp{ind} = x_s_tmp{ind}.data;
        end
        x=real(x_s_tmp{ind});
        y=imag(x_s_tmp{ind});
        x_s{ind}.x = x;
        x_s{ind}.y = y;
        
        max1=max([max(abs(x(~isinf(x)))),max(abs(y(~isinf(y)))), max1]);
      end
      zakres=[-max1 max1 -max1 max1];
    end
	end
end


% if nargin == 0, 
%   error('Not enough input arguments.'); 
% else
%   if isreal(x),
% %     y = x; x = 1:length(y); 
% %     max_y=max(abs(x(~isinf(x))));
% %     min_y=-max_x;
% %     max_x=length(y);
% %     min_x=1;
%     y = x; x = 1:length(y); 
%     max_y=10;
%     min_y=-10;
%     max_x=10;
%     min_x=10;
%   else
% %     y=imag(x);
% %     x=real(x);
% %     max_x=max(abs(x(~isinf(x))));
% %     min_x=-max_x;
% %     max_y=max_x;
% %     min_y=-max_x;
%     y=imag(x);
%     x=real(x);
%     max_x=10;
%     min_x=-10;
%     max_y=10;
%     min_y=-10;
% 
%   end
% end
if nargin==1,
  time=0.01;
else
  time=time/length(x);
end
p = 0.02;

hf=findobj(0,'tag','trajektoria');
if isempty(hf)
  hf=figure('NextPlot','add', 'Color', 'w');
else
  figure(hf);
end 
set(hf,'tag','trajektoria', 'Name', 'Trajektorie - dr in¿. Marek Blok (14.11.2006)');
delete(findobj(hf,'type', 'axes'));

global klawisz
global trajekt_hf

if pause_on_start == 0,
  klawisz = '~'; % od razu ruszaj
else
  klawisz = 13; % pauza na starcie
end
trajekt_hf = hf;
% set(hf, 'KeyPressFcn', 'global klawisz; fig = get(0, ''PointerWindow''); if fig ~= 0, klawisz = get(fig,''CurrentCharacter''); if isempty(klawisz), klawisz = ''~''; end; end;');
zoom(hf, 'off');
set(hf, 'KeyPressFcn', 'global klawisz; global trajekt_hf, fig = trajekt_hf; if fig ~= 0, klawisz = get(fig,''CurrentCharacter''); if isempty(klawisz), klawisz = ''~''; end; end;');

min_x=-5; 
max_x=5;
min_y=-5;
max_y=5;

if ile == 1,
  ax{1}=axes('Units','Normalized','Position',[0, 0, 1, 1]+ [1, 1, -2, -2]*0.05);
  if ~ishold,
    axis([min_x max_x ...
          min_y max_y]/5)
    axis square
  end
  if length(labels) >= 1
    title(labels{1}, 'Interpreter', 'latex');
  end
  set(ax{1},'Xtick', -1:0.5:1, 'Ytick', -1:0.5:1, 'Xgrid','on','Ygrid','on');
  co = get(ax{1},'colororder');
elseif ile == 2,
  for ind = 0:1,
    ax{ind+1}=axes('Units','Normalized','Position',[ind*0.5 0 [0.5 1]]+ [1, 1, -2, -2]*0.05);
    if ~ishold,
      axis([min_x max_x ...
            min_y max_y]/5)
      axis square
    end
    if length(labels) >= 2
      title(labels{ind+1}, 'Interpreter', 'latex');
    end
    set(ax{ind+1},'Xtick', -1:0.5:1, 'Ytick', -1:0.5:1, 'Xgrid','on','Ygrid','on');
  end
  co = get(ax{1},'colororder');
else
  % 12
  for ind = 0:ile-1,
    % 4 x 3
    x_ = rem(ind,4);
    y_ = (ind-x_)/4;
    ax{ind+1}=axes('Units','Normalized','Position',[x_*0.25 1-(y_+1)*0.33 [0.25 0.33]*0.95]);
    if ~ishold,
      axis([min_x max_x ...
            min_y max_y]/5)
      axis square
    end
    if length(labels) >= ile
      title(labels{ind+1}, 'Interpreter', 'latex');
    end
    set(ax{ind+1},'Xtick', -1:0.5:1, 'Ytick', -1:0.5:1, 'Xgrid','on','Ygrid','on');
  end
  co = get(ax{1},'colororder');
end



m = length(x_s{1}.x);
for ind = 1:ile-1,
  if m == 0,
    m = length(x_s{ind+1}.x);
  else
    m = min([m, length(x_s{ind+1}.x)]);
  end
end
k = ceil(p*m);

%ax1 = newplot;
% ax1=axes('Units','Pixels','Position',[30 30 300 300]);

% ax2=axes('Units','Pixels','Position',[360 30 300 300]);
for ind = 0:ile-1,
  axes(ax{ind+1})

  % Choose first three colors for head, body, and tail
  if no_head == 0,            
    head{ind+1} = line('color',co(1,:),'linestyle', 'none', 'marker', 'o', 'linewidth', 3, ...
                     'xdata',x_s{ind+1}.x(1),'ydata',x_s{ind+1}.y(1));
  end
  %body_x = []; body_y = [];
  body_x{ind+1} = NaN*ones(1,m); body_y{ind+1} = NaN*ones(1,m);
  body_x{ind+1}(1) = x_s{ind+1}.x(1);  body_y{ind+1}(1) = x_s{ind+1}.y(1);
  body{ind+1} = line('color','b','linestyle','-', 'marker','none', ... %'erase','none', ...
                   'xdata',body_x{ind+1},'ydata',body_y{ind+1});
  %tail_x = []; tail_y = [];
  tail_x{ind+1} = NaN*ones(1,k+1); tail_y{ind+1} = NaN*ones(1,k+1);
  tail_x{ind+1}(1) = x_s{ind+1}.x(1);  tail_y{ind+1}(1) = x_s{ind+1}.y(1);
  tail{ind+1} = line('color','r','linestyle','-','marker','none', ... %'erase','none', ...
                   'xdata',tail_x{ind+1},'ydata',tail_y{ind+1});
  if no_head == 0,            
    head2{ind+1} = line('color',co(4,:),'linestyle', 'none', 'marker','o', 'linewidth', 3, ...
                      'xdata',x_s{ind+1}.x(1),'ydata',x_s{ind+1}.y(1));
  end
            
  if ile > 2,
    h_punkty{ind+1} = line('color','k','linestyle', 'none', 'marker','.', 'markersize', 20, ... %'erase','none', ...
                           'xdata',0,'ydata',NaN);
  else
    h_punkty{ind+1} = line('color','k','linestyle', 'none', 'marker','.', 'markersize', 24, ... %'erase','none', ...
                           'xdata',0,'ydata',NaN);
  end;
  punkty_x{ind+1} = NaN*punkty; punkty_y{ind+1} = NaN*punkty;
  punkty_ind{ind+1} = 1;
end;

% Grow the body
for i = 2:k+1
  for ind = 0:ile-1,
    if length(x_s{ind+1}.x) > 0,
      if ~isempty(punkty)
        if any(punkty == i)
          punkty_x{ind+1}(punkty_ind{ind+1}) = x_s{ind+1}.x(i); punkty_y{ind+1}(punkty_ind{ind+1}) = x_s{ind+1}.y(i);
          punkty_ind{ind+1} = punkty_ind{ind+1} + 1;
          set(h_punkty{ind+1}, 'Xdata', punkty_x{ind+1}, 'Ydata', punkty_y{ind+1});
        end
      end
   
      j = i-1:i;
      if no_head == 0,            
        set(head{ind+1},'xdata',x_s{ind+1}.x(i),'ydata',x_s{ind+1}.y(i))
        set(head2{ind+1},'xdata',x_s{ind+1}.x(i-1),'ydata',x_s{ind+1}.y(i-1))
      end
      %set(body{ind+1},'xdata',x_s{ind+1}.x(j),'ydata',x_s{ind+1}.y(j))
      body_x{ind+1}(i) = x_s{ind+1}.x(i);  body_y{ind+1}(i) = x_s{ind+1}.y(i);
      set(body{ind+1},'xdata',body_x{ind+1},'ydata',body_y{ind+1})
      tail_x{ind+1}(i) = x_s{ind+1}.x(i);  tail_y{ind+1}(i) = x_s{ind+1}.y(i);
      set(tail{ind+1},'xdata',tail_x{ind+1},'ydata',tail_y{ind+1})
    end
  end
  
  tic;
  while toc<time
    drawnow
    pause(0)
  end;
   
  if (upper(klawisz) == 'Q') | (upper(klawisz) == 27),
    set(hf, 'KeyPressFcn', '');
    return;
  end
  if upper(klawisz) == 13,
    set(gcf, 'Color', [1, 1, 0.8]);
    klawisz = '~';
%     while length(klawisz) == 0,
    while klawisz == '~',
      drawnow
      pause(0);
    end
    set(gcf, 'Color', 'w');
    if klawisz ~= 13,
      klawisz = '~';
    end
  end
  if upper(klawisz) == '+',
    time = time / 2;
    klawisz = '~';
  end
  if upper(klawisz) == '-',
    time = time * 2;
    klawisz = '~';
  end
end


% Primary loop
for i = k+2:m
  for ind = 0:ile-1,
    if length(x_s{ind+1}.x) > 0,
      if ~isempty(punkty)
        if any(punkty == i)
          punkty_x{ind+1}(punkty_ind{ind+1}) = x_s{ind+1}.x(i); punkty_y{ind+1}(punkty_ind{ind+1}) = x_s{ind+1}.y(i);
          punkty_ind{ind+1} = punkty_ind{ind+1} + 1;
          set(h_punkty{ind+1}, 'Xdata', punkty_x{ind+1}, 'Ydata', punkty_y{ind+1});
        end
     end

     j = i-1:i;
     if no_head == 0,            
       set(head{ind+1},'xdata',x_s{ind+1}.x(i),'ydata',x_s{ind+1}.y(i))
       set(head2{ind+1},'xdata',x_s{ind+1}.x(i-1),'ydata',x_s{ind+1}.y(i-1))
     end
%      set(body{ind+1},'xdata',x_s{ind+1}.x(j),'ydata',x_s{ind+1}.y(j))
%      set(tail{ind+1},'xdata',x_s{ind+1}.x(j-k),'ydata',x_s{ind+1}.y(j-k))
     tail_x{ind+1}(1:end-1) = tail_x{ind+1}(2:end); tail_y{ind+1}(1:end-1) = tail_y{ind+1}(2:end);
     tail_x{ind+1}(end) = x_s{ind+1}.x(i);  tail_y{ind+1}(end) = x_s{ind+1}.y(i);
     set(tail{ind+1},'xdata',tail_x{ind+1},'ydata',tail_y{ind+1})
     
     body_x{ind+1}(i) = x_s{ind+1}.x(i);  body_y{ind+1}(i) = x_s{ind+1}.y(i);
     set(body{ind+1},'xdata',body_x{ind+1},'ydata',body_y{ind+1})
   end
 end
 
 tic;
 while toc<time
   drawnow
   pause(0);
 end
   
  if (upper(klawisz) == 'Q') | (upper(klawisz) == 27),
    set(hf, 'KeyPressFcn', '');
    return;
  end;
  if upper(klawisz) == 13,
    set(gcf, 'Color', [1, 1, 0.8]);
    klawisz = '~';
%     while length(klawisz) == 0,
    while klawisz == '~',
      drawnow
      pause(0);
    end;
    set(gcf, 'Color', 'w');
    if klawisz ~= 13,
      klawisz = '~';
    end
  end;
  if upper(klawisz) == '+',
    time = time / 2;
    klawisz = '~';
  end;
  if upper(klawisz) == '-',
    time = time * 2;
    klawisz = '~';
  end;
end;

% Clean up the tail
for i = m+1:m+k
  for ind = 0:ile-1,
    if length(x_s{ind+1}.x) > 0,
      j = i-1:i;
      set(tail{ind+1},'xdata',x_s{ind+1}.x(j-k),'ydata',x_s{ind+1}.y(j-k))
    end
  end;
  
  tic;
  while toc<time
   drawnow
   pause(0); %drawnow;
  end;
   
  if (upper(klawisz) == 'Q') | (upper(klawisz) == 27),
    set(hf, 'KeyPressFcn', '');
    return;
  end;
  if upper(klawisz) == 13,
    set(gcf, 'Color', [1, 1, 0.8]);
    klawisz = '~';
%     while length(klawisz) == 0,
    while klawisz == '~',
      drawnow
      pause(0);
    end;
    set(gcf, 'Color', 'w');
    if klawisz ~= 13,
      klawisz = '~';
    end
  end;
  if upper(klawisz) == '+',
    time = time / 2;
    klawisz = '~';
  end;
  if upper(klawisz) == '-',
    time = time * 2;
    klawisz = '~';
  end;
   
end;

set(hf, 'KeyPressFcn', '');

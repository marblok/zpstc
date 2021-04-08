function varargout = SPECgraf(Akcja, param)
%cos(cumsum(n/1000))

% \fixed 2017.03.22 psd ==> pwelch
% \fixed 2012.03.05 finite ==> isfinite
% \fixed 2006.10.11 dealt with LineStyle/Marker warning
% \fixed 2005.11.20 dealt with "Warning: Log of zero." in ylim evaluation
% \fixed 2005.11.20 Placed all in single file
% \fixed 2005.11.20 Fixed problems with zoom out

if nargin == 0  % LAUNCH GUI

% 	fig = openfig(mfilename,'reuse'); return; %ZapiszFig(1,'specGUI.m')
  if isempty(findobj(0, 'tag', 'Specgraf_DrawFig')) 
    
    fig = Specgraf_DrawFig;
    set(fig,'tag', 'Specgraf_DrawFig', 'Units', 'pixels');
    set(fig,'name', 'Spektrogram sygna³ów niestacjonarnych (2017.03.22) dr in¿. Marek Blok)',...
      'KeyPressFcn','1;');
    
    SPECgraf('Init');
    SPECgraf('signal');
    set(fig,'Visible','on');
%    SPECgraf('per');
  else
    SPECgraf('Exit');
  end
  return;
end;

if ~isstr(Akcja) % INVOKE NAMED SUBFUNCTION OR CALLBACK
  disp('Something''s wrong');
  return;
end;

% Generate a structure of handles to pass to callbacks, and store it. 
fig=findobj(0, 'tag', 'Specgraf_DrawFig');
hEditN=findobj(fig,'tag', 'N_edit');
hEditM=findobj(fig,'tag', 'M_edit');
hEditL=findobj(fig,'tag', 'L_edit');
hEditK=findobj(fig,'tag', 'K_edit');
hEditO=findobj(fig,'tag', 'O_edit');
hEditW=findobj(fig,'tag', 'w_edit');
hEditX=findobj(fig,'tag', 'x_edit');
h_real=findobj(fig,'tag', 'real_checkbox');
h_dB=findobj(fig,'tag', 'dB_checkbox');
hEditNoise=findobj(fig,'tag', 'Noise_edit');
hEditName=findobj(fig,'tag', 'Name_edit');
hMenu=findobj(fig,'tag', 'choose_popupmenu');
h_det=findobj(fig,'tag', 'detrend_checkbox');
ha=get(fig, 'UserData');
Ktory=get(hMenu,'Value');  
hEdit_dY=findobj(fig,'tag', 'dY_edit');

if strcmp(Akcja, 'Exit')
  close(fig);
  return;
elseif strcmp(Akcja, 'Init')
  set(hEditX,'UserData','randn(1,L)');
  set(hEditX,'String','randn(1,L)');
  set(hEditX,'Max', 2);

  set(hEditW,'UserData','boxcar(M)');
  set(hEditW,'String','boxcar(M)');

  set(hEditL,'UserData','1000');
  set(hEditL,'String','1000');
  set(hEditM,'UserData','100');
  set(hEditM,'String','100');
  set(hEditNoise,'UserData','-100');
  set(hEditNoise,'String','-100');
  
%  set(hEditN,'UserData','1');
%  set(hEditN,'String','1');
  set(hEditO,'UserData','0');
  set(hEditO,'String','0');

  set(hEditK,'UserData','256');
  set(hEditK,'String','256');

  set(hEditName,'String','new');
  set(hMenu,'String','new');

  set(h_real,'UserData',1);
  set(h_real,'Value',1);

  set(h_dB,'UserData',0);
  set(h_dB,'Value',0);
  
  set(h_det,'UserData',1);
  set(h_det,'Value',1);
  
  ha(1)=findobj(fig,'tag', 'Signal_re_axes');
  ha(2)=findobj(fig,'tag', 'Signal_im_axes');
  ha(3)=findobj(fig,'tag', 'spec_axes');
  ha(4)=findobj(fig,'tag', 'per_axes');
%  for ind=[1 2 4], axes(ha(ind)); zoom on; end;
%  axes(ha(3)); zoom off;
  
  set(fig, 'UserData', ha);

  set(hMenu,'UserData', [NaN, NaN, NaN, NaN, NaN, NaN, NaN, 1, -100, 5, -100, 5]); %handles & maximal values
  set(hEdit_dY,'String','120');
  return;
elseif strcmp(Akcja, 'change_name')
  pom=get(hMenu,'String');
  pom2=get(hEditName,'String');
  pom(Ktory,1:length(pom2)+1)=[pom2 0];
  set(hMenu,'String',pom);
  set(hMenu, 'Value', Ktory);
  return;  
elseif strcmp(Akcja, 'new')
  pom=get(hMenu,'String');
  Ktory=size(pom,1)+1;
  pom(Ktory,1:4)=['new' 0];
  set(hMenu,'String',pom);
  
  
  pom=get(hEditX,'UserData');
  pom(Ktory,1:11)=['randn(1,L)' 0];
  set(hEditX,'UserData',pom);
  pom=get(hEditW,'UserData');
  pom(Ktory,1:10)=['boxcar(M)' 0];
  set(hEditW,'UserData',pom);

  pom=get(hEditL,'UserData');
  pom(Ktory,1:4)=['100' 0];
  set(hEditL,'UserData',pom);
  pom=get(hEditM,'UserData');
  pom(Ktory,1:4)=['100' 0];
  set(hEditM,'UserData',pom);
  pom=get(hEditNoise,'UserData');
  pom(Ktory,1:5)=['-100' 0];
  set(hEditNoise,'UserData',pom);
  
%  pom=get(hEditN,'UserData');
%  pom(Ktory,1:2)=['1' 0];
%  set(hEditN,'UserData',pom);
  pom=get(hEditO,'UserData');
  pom(Ktory,1:2)=['0' 0];
  set(hEditO,'UserData',pom);

  pom=get(hEditK,'UserData');
  pom(Ktory,1:5)=['256' 0];
  set(hEditK,'UserData',pom);

  pom=get(hMenu,'UserData');
  pom(Ktory,1:size(pom,2))=ones(1,size(pom,2))*NaN;
  set(hMenu,'UserData',pom);

  pom=get(h_real,'UserData');
  pom(Ktory,1)=1;
  set(h_real,'Value', pom(Ktory,1));
  set(h_real,'UserData',pom);

  pom=get(h_det,'UserData');
  pom(Ktory,1)=1;
  set(h_det,'Value', pom(Ktory,1));
  set(h_det,'UserData',pom);

  set(hMenu,'Value', Ktory);
  
  SPECgraf('Choose');
  SPECgraf('signal', Ktory);
  return;
elseif strcmp(Akcja, 'delete')
  pom=get(hMenu,'String');
  if size(pom,1)==1,
    %SPECgraf('Reset');
    return;
  end;
  pom(Ktory,:)=[];
  set(hMenu,'String',pom);

  
  pom=get(hEditX,'UserData');
  pom(Ktory,:)=[];
  set(hEditX,'UserData',pom);
  pom=get(hEditW,'UserData');
  pom(Ktory,:)=[];
  set(hEditW,'UserData',pom);

  pom=get(hEditL,'UserData');
  pom(Ktory,:)=[];
  set(hEditL,'UserData',pom);
  pom=get(hEditM,'UserData');
  pom(Ktory,:)=[];
  set(hEditM,'UserData',pom);
  pom=get(hEditNoise,'UserData');
  pom(Ktory,:)=[];
  set(hEditNoise,'UserData',pom);
  
%  pom=get(hEditN,'UserData');
%  pom(Ktory,:)=[];
%  set(hEditN,'UserData',pom);
  pom=get(hEditO,'UserData');
  pom(Ktory,:)=[];
  set(hEditO,'UserData',pom);

  pom=get(hEditK,'UserData');
  pom(Ktory,:)=[];
  set(hEditK,'UserData',pom);

  pom=get(hMenu,'UserData');
  for ind=1:3,
    if isfinite(pom(Ktory,ind))
      delete(pom(Ktory,ind));
    end
  end
  pom(Ktory,:)=[];
  set(hMenu,'UserData',pom);

  pom=get(h_real,'UserData');
  pom(Ktory,:)=[];
  set(h_real,'UserData',pom);
  
  pom=get(h_det,'UserData');
  pom(Ktory,:)=[];
  set(h_det,'UserData',pom);
  
  set(hMenu,'Value', 1);
  
%  SPECgraf('signal', Ktory);
  SPECgraf('Choose');
  return;
  
  
elseif strcmp(Akcja, 'Choose')
  %selected new filter response
  pom=get(hMenu,'String');
  ind=find(pom(Ktory,:)==0);
  if isempty(ind)
    ind=size(pom,2);
  else
    ind=ind(1)-1;
  end
  set(hEditName,'String',pom(Ktory,1:ind));

  pom=get(hEditX,'UserData');
  set(hEditX, 'String', pom(Ktory,:));
  pom=get(hEditW,'UserData');
  set(hEditW, 'String', pom(Ktory,:));

  pom=get(hEditL,'UserData');
  set(hEditL, 'String', pom(Ktory,:));
  pom=get(hEditM,'UserData');
  set(hEditM, 'String', pom(Ktory,:));
  pom=get(hEditNoise,'UserData');
  set(hEditNoise, 'String', pom(Ktory,:));
  
%  pom=get(hEditN,'UserData');
%  set(hEditN, 'String', pom(Ktory,:));
  pom=get(hEditO,'UserData');
  set(hEditO, 'String', pom(Ktory,:));

  pom=get(hEditK,'UserData');
  set(hEditK, 'String', pom(Ktory,:));

  pom=get(h_real,'UserData');
  set(h_real,'Value', pom(Ktory,1));

  pom=get(h_det,'UserData');
  set(h_det,'Value', pom(Ktory,1));

%  SPECgraf('signal', Ktory);
  return;
  
elseif strcmp(Akcja, 'signal')
  if nargin==2,
    Ktory=param;
  else
    %save strings
    tekst=[get(hEditL, 'String') 0];
    pom=get(hEditL,'UserData');
    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
       zeros(1, max([1; size(pom,2); length(tekst)]));
    pom(Ktory,1:length(tekst))=tekst;
    set(hEditL,'UserData', setstr(pom));
    
    tekst=[get(hEditNoise, 'String') 0];
    pom=get(hEditNoise,'UserData');
    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
       zeros(1, max([1; size(pom,2); length(tekst)]));
    pom(Ktory,1:length(tekst))=tekst;
    set(hEditNoise,'UserData', setstr(pom));
    
    tekst=get(hEditX, 'String').';
    tekst=[tekst(:); 0].';
    pom=get(hEditX,'UserData');
    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
       zeros(1, max([1; size(pom,2); length(tekst)]));
    pom(Ktory,1:length(tekst))=tekst;
    set(hEditX,'UserData', setstr(pom));
    
    tekst=get(h_real, 'Value');
    pom=get(h_real,'UserData');
    pom(Ktory,1)=tekst;
    set(h_real,'UserData', pom);
  end;

  pom=get(hMenu,'UserData');
  hp_re=pom(:,1);
  hp_im=pom(:,2);
  hp_spec_t=pom(:,3);
  hp_spec_t2=pom(:,4);
  hp_spec=pom(:,5);
  hp_per=pom(:,6);
  hp_per2=pom(:,7);
  max_x=pom(:,8);
  min_spec=pom(:,9);
  max_spec=pom(:,10);
  min_per=pom(:,11);
  max_per=pom(:,12);
  
  %generate signal
  tekstL=get(hEditL,'UserData');  tekstL=tekstL(Ktory,:);
  ind=find(tekstL==0);
  if ~isempty(ind)
    tekstL=tekstL(1:ind(1)-1);
  end;
  eval(['L=' tekstL ';'], 'L=1;')

  n=0:L-1;
  tekstX=get(hEditX,'UserData');  tekstX=tekstX(Ktory,:);
  ind=find(tekstX==0);
  if ~isempty(ind)
    tekstX=tekstX(1:ind(1)-1);
  end;
  eval(['x=' tekstX ';'], 'x=zeros(1,L);')

  Re=get(h_real,'UserData');  Re=Re(Ktory,1);
  
  tekstNoise=get(hEditNoise,'UserData');  tekstNoise=tekstNoise(Ktory,:);
  ind=find(tekstNoise==0);
  if ~isempty(ind)
    tekstNoise=tekstNoise(1:ind(1)-1);
  end;
  eval(['Noise=' tekstNoise ';'], 'Noise=-300;')

  x=x(:);
  if length(x)<L;
    x(L)=0;
  else
    x=x(1:L);
  end;
  
  N_lin=10.^(Noise/20);
  if Re==1
    x=real(x)+N_lin*randn(size(x));
  else
    x=x+N_lin*(randn(size(x))+j*randn(size(x)))/sqrt(2);
  end;
  max_x(Ktory)=max(abs([real(x); imag(x)]));
   
  %draw signal

  %real part of the signal
  axes(ha(1));
  if isnan(hp_re(Ktory))
    hold on;
    hp_re(Ktory)=plot(n,real(x), 'Color', 'b');
    hold off;
  else
    set(hp_re(Ktory),'Xdata', n, 'Ydata', real(x));
  end
%  set(ha(1), 'Xlim', [-0.5, L-0.5], 'Ylim', [-1.1 1.1]*max(max_x));
%  eval('zoom reset', 'set(get(ha(1),''ZLabel''),''UserData'',[]);');    
%  reset(get(ha(1),'ZLabel'));    
  
  %imaginary part of the signal
%  axes(ha(2));
  if isnan(hp_im(Ktory))
    hold on;
    hp_im(Ktory)=plot(n,imag(x), 'Color', 'r');
    hold off;
  else
    set(hp_im(Ktory),'Xdata', n, 'Ydata', imag(x));
  end
  set(ha(1), 'Xlim', [-0.5, L-0.5], 'Ylim', [-1.1 1.1]*max(max_x));
%  set(get(ha(2),'ZLabel'),'UserData',[]);    
%  reset(get(ha(2),'ZLabel'));    
%  eval('zoom reset', 'set(get(ha(1),''ZLabel''),''UserData'',[]);');    
  if L>512
    set([hp_re, hp_im], 'Marker', '.', 'MarkerSize', 4);
  else
    set([hp_re, hp_im], 'LineStyle', '-');
  end;
  eval('rmappdata(get(ha(1),''Zlabel''),''ZOOMAxesData'')','set(get(ha(1),''ZLabel''),''UserData'',[])');

  set(hMenu,'UserData', [hp_re, hp_im,  hp_spec_t, hp_spec_t2, hp_spec, hp_per, hp_per2, max_x, min_spec, max_spec, min_per, max_per]);

  %compute and draw periodogram
  SPECgraf('spec', Ktory)
%  SPECgraf zoom_on;
  return;
  
elseif strcmp(Akcja, 'spec')
  if nargin==2,
    Ktory=param;
  else
    %save strings
    tekst=[get(hEditK, 'String') 0];
    pom=get(hEditK,'UserData');
    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
       zeros(1, max([1; size(pom,2); length(tekst)]));
    pom(Ktory,1:length(tekst))=tekst;
    set(hEditK,'UserData', setstr(pom));
    
    tekst=[get(hEditM, 'String') 0];
    pom=get(hEditM,'UserData');
    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
       zeros(1, max([1; size(pom,2); length(tekst)]));
    pom(Ktory,1:length(tekst))=tekst;
    set(hEditM,'UserData', setstr(pom));
    
%    tekst=[get(hEditN, 'String') 0];
%    pom=get(hEditN,'UserData');
%    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
%       zeros(1, max([1; size(pom,2); length(tekst)]));
%    pom(Ktory,1:length(tekst))=tekst;
%    set(hEditN,'UserData', setstr(pom));
    
    tekst=[get(hEditO, 'String') 0];
    pom=get(hEditO,'UserData');
    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
       zeros(1, max([1; size(pom,2); length(tekst)]));
    pom(Ktory,1:length(tekst))=tekst;
    set(hEditO,'UserData', setstr(pom));
    
    tekst=[get(hEditW, 'String') 0];
    pom=get(hEditW,'UserData');
    pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
       zeros(1, max([1; size(pom,2); length(tekst)]));
    pom(Ktory,1:length(tekst))=tekst;
    set(hEditW,'UserData', setstr(pom));
    
    tekst=get(h_det, 'Value');
    pom=get(h_det,'UserData');
    pom(Ktory,1)=tekst;
    set(h_det,'UserData', pom);
  end;

  pom=get(hMenu,'UserData');
  hp_re=pom(:,1);
  hp_im=pom(:,2);
  hp_spec_t=pom(:,3);
  hp_spec_t2=pom(:,4);
  hp_spec=pom(:,5);
  hp_per=pom(:,6);
  hp_per2=pom(:,7);
  max_x=pom(:,8);
  min_spec=pom(:,9);
  max_spec=pom(:,10);
  min_per=pom(:,11);
  max_per=pom(:,12);
  
  %generate signal
  tekstK=get(hEditK,'UserData');  tekstK=tekstK(Ktory,:);
  ind=find(tekstK==0);
  if ~isempty(ind)
    tekstK=tekstK(1:ind(1)-1);
  end;
  eval(['K=' tekstK ';'], 'K=16;')

  tekstM=get(hEditM,'UserData');  tekstM=tekstM(Ktory,:);
  ind=find(tekstM==0);
  if ~isempty(ind)
    tekstM=tekstM(1:ind(1)-1);
  end;
  eval(['M=' tekstM ';'], 'M=16;')
  if M>K
    M=K;
    set(hEditM,'String', num2str(M));
  end;

%  tekstN=get(hEditN,'UserData');  tekstN=tekstN(Ktory,:);
%  ind=find(tekstN==0);
%  if ~isempty(ind)
%    tekstN=tekstN(1:ind(1)-1);
%  end;
%  eval(['N=' tekstN ';'], 'N=16;')

  tekstO=get(hEditO,'UserData');  tekstO=tekstO(Ktory,:);
  ind=find(tekstO==0);
  if ~isempty(ind)
    tekstO=tekstO(1:ind(1)-1);
  end;
%   eval(['O=' tekstO ';'], 'O=0;')
  O=eval(tekstO, '0'); % \Fixed 2005.11.03 
  O=round(O/100*M); % \Fixed 2005.11.03 nak³adkowanie podawane w procentach !!! 
  

  tekstW=get(hEditW,'UserData');  tekstW=tekstW(Ktory,:);
  ind=find(tekstW==0);
  if ~isempty(ind)
    tekstW=tekstW(1:ind(1)-1);
  end;
  eval(['w=' tekstW ';'], 'w=ones(1,M);')

  w=w(:);
  if length(w)<M;
    w(M)=0;
  else
    w=w(1:M);
  end;
  
  x=get(hp_re(Ktory), 'Ydata')+j*get(hp_im(Ktory), 'Ydata');

%   O=floor(O/100*M); % \Fixed 2005.11.03
  if O>=M
    O=M-1;
    set(hEditO,'String', num2str(O/M*100));
  end;
  N = fix((length(x)-O)/(M-O));
  set(hEditN, 'String', sprintf('%i', N));

%  det_=get(h_det,'UserData');  det_=det_(Ktory,1);
%  if det_==0,
%    det='none';
%  else
%    det_='linear';
%  end;
  
  [Spec, f, t]=specgram(x, K, 1, w, O);
  Spec=(abs(Spec).^2)/norm(w)^2;
%   [Per, f2]=psd(x, K, 1, w, O);
  [Per, f2]=pwelch(x, w, O, K, 1);
  Spec=abs(Spec);

  Re=get(h_real,'UserData'); 
  if ~Re(Ktory,1)
    f=fftshift(f);
    f(1:floor(K/2))=f(1:floor(K/2))-1;
    f2=fftshift(f2);
    f2(1:floor(K/2))=f2(1:floor(K/2))-1;
    Spec=[Spec(ceil(K/2):K,:); Spec(1:floor(K/2),:)];
    Per=fftshift(Per);
  end
  
  min_spec(Ktory)=min([min(Spec(isfinite(Spec))); 0]);
  max_spec(Ktory)=max([max(Spec(isfinite(Spec))); 0.001]);
  min_per(Ktory)=min([min(Per(isfinite(Per))); 0]);
  max_per(Ktory)=max([max(Per(isfinite(Per))); 0.001]);
  if get(h_dB, 'UserData') %dB
    Spec=10*log10(Spec);
    Per=10*log10(Per);
  end
  set(get(ha(2),'Ylabel'),'UserData', f);
  set(get(ha(3),'Ylabel'),'UserData', Spec);
  set(get(ha(4),'Ylabel'),'UserData', Per);
  
  %draw signal
  axes(ha(2));
  if length(t)>1
    t2=t+(t(2)-t(1))/2;
  else
    t2=M/2;
  end;
  
  if isnan(hp_spec_t(Ktory))
    hold on;
%    hp_spec_t(Ktory)=plot(t2,max(Spec), 'Color', 'k');
    pomoc=abs(get(hp_re(Ktory),'Ydata')+j*get(hp_im(Ktory),'Ydata'));
    hp_spec_t(Ktory)=plot(get(hp_re(Ktory),'Xdata'), pomoc, 'Color', 'k');
    if length(pomoc)>512
     set(hp_spec_t(Ktory), 'Marker', '.', 'MarkerSize', 4);
   else
     set(hp_spec_t(Ktory), 'LineStyle', '-');
   end;
    hold off;
  else
%    set(hp_spec_t(Ktory),'Xdata', t2, 'Ydata', max(Spec), 'Color', 'k');
    set(hp_spec_t(Ktory),'Xdata', get(hp_re(Ktory),'Xdata'),...
                         'Ydata', abs(get(hp_re(Ktory),'Ydata')+j*get(hp_im(Ktory),'Ydata')), 'Color', 'k');
  end
  if length(t)==1,
    set(ha(2),  'Xlim', [t(1)-0.5 t(1)+0.5], 'Ylim', [-1.1 1.1]*max(max_x));
  else
    set(ha(2),  'Xlim', [t(1) t(length(t))+(t(2)-t(1))], 'Ylim', [-1.1 1.1]*max(max_x));
  end;
%  set(get(ha(2),'ZLabel'),'UserData',[]);    
%  reset(get(ha(2),'ZLabel'));    
  eval('rmappdata(get(ha(2),''Zlabel''),''ZOOMAxesData'')','set(get(ha(2),''ZLabel''),''UserData'',[])');

  %spektrogram
  axes(ha(3));
  if isnan(hp_spec(Ktory))
    hold on;
    hp_spec(Ktory)=image(t2, f, Spec);
    colormap(hot);
    hold off;
  else
    set(hp_spec(Ktory),'Xdata', t2, 'Ydata', f, 'Cdata', Spec);
  end
  if length(t)==1,
    tlim_=[t(1)-0.5 t(1)+0.5];
  else
    tlim_=[t(1) t(length(t))+t(2)-t(1)];
  end;
  if all(Re)
    set(ha(3), 'Ylim', [0 0.5], 'Xlim', tlim_);
  else
    set(ha(3), 'Ylim', [-0.5 0.5], 'Xlim', tlim_);
  end
%  set(get(ha(3),'ZLabel'),'UserData',[]);    
%  reset(get(ha(3),'ZLabel'));    
  eval('zoom reset', 'set(get(ha(3),''ZLabel''),''UserData'',[]);');    
  
  axes(ha(4));
  if isnan(hp_per(Ktory))
    hold on;
    hp_per(Ktory)=plot(Per, f2, 'Color', 'k');
    hold off;
  else
    set(hp_per(Ktory),'Ydata', f2, 'Xdata', Per, 'Color', 'k');
  end
  set(ha(4), 'Xdir', 'reverse', 'YTick', []); %, 'Xlim', [t(1) t(length(t))], 'Ylim', [-1.1 1.1]*max(max_x));
  if all(Re)
    set(ha(4), 'Ylim', [0 0.5]);
  else
    set(ha(4), 'Ylim', [-0.5 0.5]);
  end
  eval('rmappdata(get(ha(4),''Zlabel''),''ZOOMAxesData'')','set(get(ha(4),''ZLabel''),''UserData'',[])');

  set(hMenu,'UserData', [hp_re, hp_im, hp_spec_t, hp_spec_t2, hp_spec, hp_per, hp_per2, max_x, min_spec, max_spec, min_per, max_per]);

  %  delete
  if 1
    pom=get(hMenu,'UserData');
    hp_=pom(:,[4,7]);
    hp_=hp_(isfinite(hp_));
    if length(hp_)>0
      delete(hp_)
      pom(:,4)=NaN;
      pom(:,7)=NaN;
      set(hMenu,'UserData',pom);
    end;
  end
  
  SPECgraf('spec_ylim');
  SPECgraf zoom_on;
  SPECgraf zoom_off;
  return;
  
elseif strcmp(Akcja, 'dB')
  pom=get(h_dB, 'UserData');
  if pom~=get(h_dB, 'Value')
    pom=get(hMenu,'UserData');
    hp_spec=pom(:,3);
    for ind=1:length(hp_spec)
      sygn=get(get(ha(3),'Ylabel'),'UserData');
      per=get(get(ha(4),'Ylabel'),'UserData');
      if get(h_dB, 'Value') %dB
        sygn=10*log10(sygn);
        per=10*log10(per);
      else %lin
        sygn=10.^(sygn/10);
        per=10.^(per/10);
      end;
      set(get(ha(3),'Ylabel'),'UserData', sygn);
      set(get(ha(4),'Ylabel'),'UserData', per);
%      set(hp_spec(ind),'Cdata', sygn);
    end;
    set(h_dB, 'UserData', get(h_dB, 'Value'));

    hp_=pom(:,[4 7]);
    hp_=hp_(isfinite(hp_));
    if length(hp_)>0
      delete(hp_)
      pom(:,4)=NaN;
      pom(:,7)=NaN;
      set(hMenu,'UserData',pom);
    end;
    SPECgraf('spec_ylim');
  end
  return
  
elseif strcmp(Akcja, 'spec_ylim')
  pom=get(hMenu,'UserData');
  hp_re=pom(:,1);
  hp_im=pom(:,2);
  hp_spec_t=pom(:,3);
  hp_spec=pom(:,5);
  hp_per=pom(:,6);
  min_spec=pom(:,9);
  max_spec=pom(:,10);
  min_per=pom(:,11);
  max_per=pom(:,12);
  if get(h_dB, 'UserData') %dB
    tekst=get(hEdit_dY,'String');
    eval(['dY=' tekst ';'], 'dY=120;');
    if dY<=0, dY=10; end;
    params_=[min(min_spec) max(max_spec)];
    ind_params = find(abs(params_) <= eps);
    if length(ind_params) > 0,
      params_(ind_params) = NaN*ones(size(ind_params));
    end
    ylim_=10*log10(params_);
    if ~isfinite(ylim_(2))
      ylim_(2)=0;
    end
    ylim_(1)=ylim_(2)-dY;
    params_=[min(min_per) max(max_per)];
    ind_params = find(abs(params_) <= eps);
    if length(ind_params) > 0,
      params_(ind_params) = NaN*ones(size(ind_params));
    end
    ylim_per=10*log10(params_);
    if ~isfinite(ylim_per(2))
      ylim_per(2)=0;
    end
    ylim_per(1)=ylim_per(2)-dY;
  else
    ylim_=[0 max(max_spec)];
    ylim_per=[0 max(max_per)];
  end
  ylim_per(2)=max([ylim_per(2) ylim_(2)]);
  f=get(get(ha(2),'Ylabel'),'UserData');
  Spec=get(get(ha(3),'Ylabel'),'UserData');
  Per=get(get(ha(4),'Ylabel'),'UserData');
  set(hp_spec_t(Ktory),'Ydata', abs(get(hp_re(Ktory),'Ydata')+j*get(hp_im(Ktory),'Ydata')));
  set(ha(2),'Ylim', ylim_+(ylim_(2)-ylim_(1))*[-0.1 0.1]);
  Spec=64*(Spec-ylim_(1))/(ylim_(2)-ylim_(1));

  set(hp_spec(Ktory),'Cdata', Spec);
  set(hp_per(Ktory),'Xdata', Per, 'Ydata', f);
  set(ha(4),'Xlim', ylim_per);

  if get(h_dB, 'UserData') %dB
    set(hp_spec_t(Ktory),'Ydata', 20*log10(abs(get(hp_re(Ktory),'Ydata')+j*get(hp_im(Ktory),'Ydata'))));
  else
    set(hp_spec_t(Ktory),'Ydata', abs(get(hp_re(Ktory),'Ydata')+j*get(hp_im(Ktory),'Ydata')));
  end
  
  SPECgraf zoom_on;
  eval('zoom reset', 'set(get(ha(3),''ZLabel''),''UserData'',[]);');    
  SPECgraf zoom_off;
  return
elseif strcmp(Akcja, 'zoom_on')
  zoom on;
%   pom=get(fig, 'WindowButtonDownFcn');
%   set(fig, 'WindowButtonDownFcn', 'SPECgraf zoom');
%   set(get(ha(3),'Xlabel'), 'UserData', pom);
  return;  
elseif strcmp(Akcja, 'zoom_off')
  if get(findobj(fig,'tag','checkbox_zoom'),'Value') == 0,
%     pom = get(get(ha(3),'Xlabel'), 'UserData');
%     set(fig, 'WindowButtonDownFcn', pom);
    zoom off;
    set(fig, 'WindowButtonDownFcn', 'SPECgraf zoom');
    set(get(ha(3),'Xlabel'), 'UserData', '1;');
  end
  return;  
elseif strcmp(Akcja, 'zoom_spec')
  if get(findobj(fig,'tag','checkbox_zoom'),'Value') ~= 0,
    Specgraf zoom_on;
  else
    Specgraf zoom_off;
  end
elseif strcmp(Akcja, 'zoom')
%  if strcmp(get(fig,'SelectionType'),'normal') | (gca~=ha(3))
  pause(0);
  if (get(gco,'Parent')~=ha(3)) | get(findobj(fig,'tag','checkbox_zoom'),'Value')
    eval(get(get(ha(3),'Xlabel'), 'UserData'));
  elseif get(gco,'Parent')==ha(3)
    pom=get(ha(3), 'CurrentPoint');
    f_=pom(1,2); t_=pom(1,1);

    pom_=get(hMenu,'UserData');
    hp_spec_t2=pom_(:,4);
    hp_spec=pom_(:,5);
    hp_per2=pom_(:,7);
    f=get(hp_spec(Ktory), 'Ydata');
    ind_f=find(abs(f-f_)==min(abs(f-f_))); ind_f=ind_f(1);
    t=get(hp_spec(Ktory), 'Xdata');
%    if length(t)>1,
%      t=t+(t(2)-t(1))/2;
%    end;
    ind_t=find(abs(t-t_)==min(abs(t-t_)));
    ind_t=ind_t(1);
    set(findobj(fig,'tag', 'text_t_f'),...
      'ForegroundColor', 'r', 'String', sprintf('n=%i, f=%6.3f', t(ind_t),f(ind_f)));

    Spec=get(get(ha(3),'Ylabel'),'UserData');
    axes(ha(4));
    if length(Spec(:,ind_t))>length(f)
      ind=find(f==0);
      Spec(ind(1),:)=[];
    end;
    if isnan(hp_per2(Ktory))
      hold on;
      hp_per2(Ktory)=plot(Spec(:,ind_t),f, 'Color', 'r', 'LineWidth', 2);
      hold off;
    else
      set(hp_per2(Ktory),'Xdata', Spec(:,ind_t), 'Ydata', f);
    end
%    pom=get(ha(4), 'Xlim');
%    pom(2)=max([pom(2) max(Spec(:,ind_t))]);
%    set(ha(4),'Xlim', pom);

    axes(ha(2));
    if isnan(hp_spec_t2(Ktory))
      hold on;
      hp_spec_t2(Ktory)=plot(t,Spec(ind_f,:), 'Color', 'r', 'LineWidth', 2);
      hold off;
      SPECgraf zoom_on;
      SPECgraf zoom_off;
    else
      set(hp_spec_t2(Ktory),'Xdata', t, 'Ydata', Spec(ind_f,:));
    end

    pom_(Ktory,7)=hp_per2;
    pom_(Ktory,4)=hp_spec_t2;
    set(hMenu,'UserData',pom_);
  end;
 return;  
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig_h=Specgraf_DrawFig
fig_h=figure;
set(fig_h, ...
  'Units', 'Normalized',...
  'Position', [0.1125 0.1188 0.7813 0.7448],...
  'Color', [1.0000 1.0000 1.0000],...
  'Name', 'Untitled',...
  'NumberTitle', 'off',...
  'Menu', 'none',...
  'Tag', 'figure1'...
  );
h=axes;
set(h, ...
  'Units', 'Normalized',...
  'Position', [0.0330 0.0881 0.2150 0.4867],...
  'Color', [1.0000 1.0000 1.0000],...
  'XColor', [0.0000 0.0000 0.0000],...
  'YColor', [0.0000 0.0000 0.0000],...
  'FontSize', 10,...
  'Box', 'on',...
  'Tag', 'per_axes'...
  );
h=axes;
set(h, ...
  'Units', 'Normalized',...
  'Position', [0.2900 0.6350 0.6990 0.1538],...
  'Color', [1.0000 1.0000 1.0000],...
  'XColor', [0.0000 0.0000 0.0000],...
  'YColor', [0.0000 0.0000 0.0000],...
  'FontSize', 10,...
  'Box', 'on',...
  'Tag', 'Signal_im_axes'...
  );
h=axes;
set(h, ...
  'Units', 'Normalized',...
  'Position', [0.2880 0.0867 0.7030 0.4867],...
  'Color', [1.0000 1.0000 1.0000],...
  'XColor', [0.0000 0.0000 0.0000],...
  'YColor', [0.0000 0.0000 0.0000],...
  'FontSize', 10,...
  'Box', 'on',...
  'Tag', 'spec_axes'...
  );
h=axes;
set(h, ...
  'Units', 'Normalized',...
  'Position', [0.2920 0.8266 0.6980 0.1580],...
  'Color', [1.0000 1.0000 1.0000],...
  'XColor', [0.0000 0.0000 0.0000],...
  'YColor', [0.0000 0.0000 0.0000],...
  'FontSize', 10,...
  'Box', 'on',...
  'Tag', 'Signal_re_axes'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'frame',...
  'Units', 'Normalized',...
  'Position', [0.0050 0.7972 0.2400 0.1930],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', '',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'frame2'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'frame',...
  'Units', 'Normalized',...
  'Position', [0.0050 0.5916 0.2400 0.2028],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', '',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'frame4'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.1190 0.6364 0.0280 0.0350],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', '[%]',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text23'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'checkbox',...
  'Units', 'Normalized',...
  'Position', [0.8720 0.0126 0.0660 0.0196],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'zoom',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'Specgraf zoom_spec',...
  'Tag', 'checkbox_zoom'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'pushbutton',...
  'Units', 'Normalized',...
  'Position', [0.0140 0.0098 0.0660 0.0378],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'Exit',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf Exit',...
  'Tag', 'pushbutton5'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.4640 0.0028 0.3540 0.0350],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', '',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text_t_f'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.1750 0.0126 0.0750 0.0364],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''spec_ylim'')',...
  'Tag', 'dY_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.1110 0.0154 0.0620 0.0308],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'dY [dB] =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text21'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'checkbox',...
  'Units', 'Normalized',...
  'Position', [0.9270 0.5748 0.0500 0.0322],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'dB',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''dB'')',...
  'Tag', 'dB_checkbox'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'pushbutton',...
  'Units', 'Normalized',...
  'Position', [0.0660 0.4909 0.0600 0.0420],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'Delete',...
  'Value', 0,...
  'Visible', 'off',...
  'Callback', 'SPECgraf(''delete'')',...
  'Tag', 'pushbutton4'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'pushbutton',...
  'Units', 'Normalized',...
  'Position', [0.0120 0.4909 0.0480 0.0420],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'New',...
  'Value', 0,...
  'Visible', 'off',...
  'Callback', 'SPECgraf(''new'')',...
  'Tag', 'pushbutton3'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.1580 0.6769 0.0750 0.0364],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''spec'')',...
  'Tag', 'N_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.1280 0.6825 0.0300 0.0266],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'N =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text18'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.0390 0.8014 0.0720 0.0378],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''signal'')',...
  'Tag', 'L_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0090 0.8028 0.0290 0.0266],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'L =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text17'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'checkbox',...
  'Units', 'Normalized',...
  'Position', [0.1290 0.7301 0.1100 0.0252],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'detrend',...
  'Value', 0,...
  'Visible', 'off',...
  'Callback', 'SPECgraf(''spec'')',...
  'Tag', 'detrend_checkbox'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.0400 0.6378 0.0750 0.0364],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''spec'')',...
  'Tag', 'O_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0110 0.6378 0.0280 0.0294],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'O =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text16'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.0520 0.6000 0.1870 0.0336],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''spec'')',...
  'Tag', 'w_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0130 0.5986 0.0400 0.0294],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'w[n] =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text15'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.1670 0.8042 0.0720 0.0364],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''signal'')',...
  'Tag', 'Noise_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.1190 0.8042 0.0470 0.0266],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'Noise=',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text14'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.0410 0.6769 0.0750 0.0364],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''spec'')',...
  'Tag', 'M_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0110 0.6825 0.0300 0.0266],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'M =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text13'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.1340 0.4909 0.1030 0.0392],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'off',...
  'Callback', 'SPECgraf(''change_name'')',...
  'Tag', 'Name_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'popupmenu',...
  'Units', 'Normalized',...
  'Position', [0.0120 0.5413 0.2260 0.0294],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Popup Menu',...
  'Value', 1,...
  'Visible', 'off',...
  'Callback', 'SPECgraf(''Choose'')',...
  'Tag', 'choose_popupmenu'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'checkbox',...
  'Units', 'Normalized',...
  'Position', [0.1160 0.9580 0.1110 0.0266],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'real signal',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''signal'')',...
  'Tag', 'real_checkbox'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0100 0.7622 0.1470 0.0266],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'Spectrograf parameters',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text10'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.0400 0.7259 0.0750 0.0364],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''spec'')',...
  'Tag', 'K_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0140 0.7329 0.0240 0.0252],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'K =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text9'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'edit',...
  'Units', 'Normalized',...
  'Position', [0.0390 0.8462 0.2010 0.1077],...
  'BackGroundColor', [1.0000 1.0000 1.0000],...
  'String', 'Edit Text',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', 'SPECgraf(''signal'')',...
  'Tag', 'x_edit'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0130 0.9245 0.0240 0.0252],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'x =',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text4'...
  );
h=uicontrol;
x=version;
if (x(1)>='5'),
  set(h, ...
    'FontSize', 10);
end;
set(h, ...
  'Style', 'text',...
  'Units', 'Normalized',...
  'Position', [0.0120 0.9636 0.1010 0.0238],...
  'BackGroundColor', [0.9255 0.9137 0.8471],...
  'String', 'Testing signal',...
  'Value', 0,...
  'Visible', 'on',...
  'Callback', '',...
  'Tag', 'text2'...
  );

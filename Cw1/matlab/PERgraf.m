function PERgraf(Akcja, param)

global supress_msgbox

% \fixed 2018.03.02 pwelch: warning for L < M
% \fixed 2015.10.19 pwelch + detrend is now used instead of psd
%
% \fixed 2006.12.20 segment length for raw and modified periodographs is
%                   now automaticaly set based on signal length
%                   This spares us annoying warnings about L and M differences
% \added 2006.10.10 added 'lock' option to lock zoom scale
% \added 2006.10.09 added Fs, it can also be used in edit fields like 'x' in the same way like A 
% \added 2006.10.09 added popupmenu with periodogrph type selection
% \fixed 2005.11.10 "delete" for last element resets it
% \fixed 2005.11.10 "new" now copies current setting instead off standard ones 
% \Fixed 2005.11.10 dodano "supress_msgbox"
% \Fixed 2005.11.10 przy zmianie A wszystkie wykresy s¹ przerysowane
% \Fixed 2005.11.10 niewidoczny wspó³czynnik 10^-x dla ma³ych sygna³ów
% \Fixed 2005.11.10 dodano akcjê "actionA"
%   - "A_edit"
%   - "dA_edit"
%   - "A_slider"
%   - "minA_edit"
%   - "maxA_edit"
%   - "A_plus_radio"
%   - "A_razy_radio"
%
% \Fixed 2005.11.02 nak³adkowanie podawane w procentach !!! 
% \Fixed 2005.11.02 detrender domyœlnie wy³¹czony
% \Fixed 2005.11.02 Prze³¹czanie pomiêdzy wzorcami nie uaktualnia N
% \Fixed 2005.11.02 message when error in entered data
% \Fixed 2005.11.02 wyœwietlany fragment sygna³u o d³ugoœci ostanio
%   wprowadzonego L, a nie najd³u¿szego z wszystkich sygna³ów (added max_L)
% \Fixed 2005.11.02 Refresh Button
% \Fixed 2005.11.02 Use hourglass when busy
% \Fixed 2005.11.02 Mean, var should be updated acording to choose and use colour
% \Fixed 2005.11.02 Warning when M > K.
%
%
% Szacowanie d³ugoœci segmentu na podstawie liczby segmentów N, d³ugoœci ci¹gu L i wielkoœci nak³adkowania O: M = L/(N-O*(N-1))
if nargin == 0  % LAUNCH GUI

%  fig = openfig(mfilename,'reuse'); return; %ZapiszFig(1,'PerGUI.m')
  if isempty(findobj(0, 'tag', 'perGUI')) 
%    set(0, 'DefaultFigureVisible', 'off'); 	   
    fig = perGUI;
%    set(0, 'DefaultFigureVisible', 'on'); 	   
    set(fig,'tag', 'perGUI', 'Visible', 'on', 'Units', 'pixels');
    set(fig,'NumberTitle', 'off', 'Name','Periodogramowe estymatory widma gêstoœci mocy ver. 1.3e (dr in¿. Marek Blok 02.03.2018)')    

    PERgraf('Init');
    PERgraf('signal');
%    PERgraf('per');
  else
    PERgraf('Exit');
  end
  return;
end;

if ~isstr(Akcja) % INVOKE NAMED SUBFUNCTION OR CALLBACK
  disp('Something''s wrong');
  return;
end;

% Generate a structure of handles to pass to callbacks, and store it. 
fig=findobj(0, 'tag', 'perGUI');
hEditN=findobj(fig,'tag', 'N_edit');
hEditM=findobj(fig,'tag', 'M_edit');
hEditL=findobj(fig,'tag', 'L_edit');
hEditK=findobj(fig,'tag', 'K_edit');
hEditO=findobj(fig,'tag', 'O_edit');
hEditW=findobj(fig,'tag', 'w_edit');
hEditX=findobj(fig,'tag', 'x_edit');
h_real=findobj(fig,'tag', 'real_checkbox');
h_dB=findobj(fig,'tag', 'dB_checkbox');
h_rgb(1)=findobj(fig,'tag', 'r_checkbox');
h_rgb(2)=findobj(fig,'tag', 'g_checkbox');
h_rgb(3)=findobj(fig,'tag', 'b_checkbox');
hEditNoise=findobj(fig,'tag', 'Noise_edit');
hEditName=findobj(fig,'tag', 'Name_edit');
hMenu=findobj(fig,'tag', 'choose_popupmenu');
h_det=findobj(fig,'tag', 'detrend_checkbox');
ha=get(fig, 'UserData');
Ktory=get(hMenu,'Value');  
hEdit_dY=findobj(fig,'tag', 'dY_edit');
hParam=findobj(fig,'tag', 'text20');
% \fixed 2005.11.10
hA.A_edit=findobj(fig,'tag', 'A_edit');
hA.dA_edit=findobj(fig,'tag', 'dA_edit');
hA.A_slider=findobj(fig,'tag', 'A_slider');
hA.minA_edit=findobj(fig,'tag', 'minA_edit');
hA.maxA_edit=findobj(fig,'tag', 'maxA_edit');
hA.A_plus_radio=findobj(fig,'tag', 'A_plus_radio');
hA.A_razy_radio=findobj(fig,'tag', 'A_razy_radio');
% /fixed 2005.11.10
% \added 2006.10.09
h_per_type=findobj(fig,'tag', 'per_type_popupmenu');
h_Fs=findobj(fig,'tag', 'Fs_edit');
h_lock=findobj(fig,'tag', 'lock_checkbox');
% /added 2006.10.09


if strcmp(Akcja, 'Exit')
  close(fig);
  return;
elseif strcmp(Akcja, 'Init')
  % UserData przechowuje wartoœci domyœlne
  supress_msgbox = 0;
  
  set(hEditX,'UserData','randn(1,L)');
  set(hEditX,'String','randn(1,L)');
  set(hEditX,'Max', 2);

  set(hEditW,'UserData','boxcar(M)');
  set(hEditW,'String','boxcar(M)');

  set(hEditL,'UserData','100');
  set(hEditL,'String','100');
  set(hEditM,'UserData','100');
  set(hEditM,'String','100');
  set(hEditNoise,'UserData','-600');
  set(hEditNoise,'String','-600');
  
%  set(hEditN,'UserData','1');
%  set(hEditN,'String','1');
  set(hEditO,'UserData','0');
  set(hEditO,'String','0');

  set(hEditK,'UserData','1024');
  set(hEditK,'String','1024');

  set(hEditName,'String','new');
  set(hMenu,'String','new');

  set(h_real,'UserData',1);
  set(h_real,'Value',1);

  set(h_dB,'UserData',0);
  set(h_dB,'Value',0);
  
  set(h_det,'UserData',0);  % \Fixed 2005.11.02 - detrender domyœlnie wy³¹czony
  set(h_det,'Value',0);
  
  pom(1,1:3)=[0 0 0];
  for ind=1:3,
    set(h_rgb(ind),'Value', pom(1,ind));
  end;
  set(h_rgb(1),'UserData',pom);

  ha(1)=findobj(fig,'tag', 'Signal_re_axes');
  ha(2)=findobj(fig,'tag', 'Signal_im_axes');
  ha(3)=findobj(fig,'tag', 'per_axes');
  set(fig, 'UserData', ha, 'Units', 'Pixel');

  set(hMenu,'UserData', [NaN, NaN, NaN, 1, -100, 5, 1, 0, 0]); %handles & maximal values

  set(hEdit_dY,'String','120');
  
  % \fixed 2005.11.10
  A.A = 1.0; A.dA = 0.1;
  A.min=0.0; A.max=2.0;
  A.is_plus=1;
  set(hA.A_edit, 'UserData', A);
  % UpdateA(hA);
  set(hA.A_edit, 'String', A.A);
  set(hA.dA_edit, 'String', A.dA);
  set(hA.minA_edit, 'String', A.min);
  set(hA.maxA_edit, 'String', A.max);
  set(hA.A_plus_radio, 'Value', A.is_plus);
  set(hA.A_razy_radio, 'Value', 1-A.is_plus);
  set(hA.A_slider, 'min', A.min, 'max', A.max, 'value', A.A);
  % /fixed 2005.11.10
  % \added 2006.10.09
  set(h_per_type, 'UserData', 1); % surowy periodograf jako domyœlny
  set(h_Fs, 'UserData', 'Inf');
  set(h_Fs, 'String', 'Inf');
  % \added 2006.10.09
  
  UpdatePergraphGUI(1, fig);
  return;
  
elseif strcmp(Akcja, 'change_name')
  pom=get(hMenu,'String');
  pom2=get(hEditName,'String');
  pom(Ktory,1:length(pom2)+1)=[pom2 0];
  set(hMenu,'String',pom);
  set(hMenu, 'Value', Ktory);
  return;  

elseif strcmp(Akcja, 'Fs_changed')
  tekst=[get(h_Fs, 'String') 0];
  pom=get(h_Fs,'UserData');
  pom(Ktory,1:max([1; size(pom,2); length(tekst)]))= ...
      zeros(1, max([1; size(pom,2); length(tekst)]));
  pom(Ktory,1:length(tekst))=tekst;
  set(h_Fs,'UserData', setstr(pom));
  
  tekstX=get(hEditX,'UserData');  tekstX=tekstX(Ktory,:);
  ind=find(tekstX==0);
  if ~isempty(ind)
    tekstX=tekstX(1:ind(1)-1);
  end;
  if ~isempty(findstr(tekstX, 'Fs'))
    PERgraf('signal', Ktory);
  end
  PERgraf('per', Ktory);

elseif strcmp(Akcja, 'new')
  old_Ktory=get(hMenu,'Value');  

  pom=get(hMenu,'String');
  Ktory=size(pom,1)+1;
  pom(Ktory,1:4)=['new' 0];
  set(hMenu,'String',pom);
  
  
  pom=get(hEditX,'UserData');
%   pom(Ktory,1:11)=['randn(1,L)' 0];
  pom(Ktory,:)=pom(old_Ktory,:);
  set(hEditX,'UserData',pom);
  pom=get(hEditW,'UserData');
%   pom(Ktory,1:10)=['boxcar(M)' 0];
  pom(Ktory,:)=pom(old_Ktory,:);
  set(hEditW,'UserData',pom);

  pom=get(hEditL,'UserData');
%   pom(Ktory,1:4)=['100' 0];
  pom(Ktory,:)=pom(old_Ktory,:);
  set(hEditL,'UserData',pom);
  pom=get(hEditM,'UserData');
%   pom(Ktory,1:4)=['100' 0];
  pom(Ktory,:)=pom(old_Ktory,:);
  set(hEditM,'UserData',pom);
  pom=get(hEditNoise,'UserData');
%   pom(Ktory,1:5)=['-100' 0];
  pom(Ktory,:)=pom(old_Ktory,:);
  set(hEditNoise,'UserData',pom);
  
%  pom=get(hEditN,'UserData');
%  pom(Ktory,1:2)=['1' 0];
%  set(hEditN,'UserData',pom);
  pom=get(hEditO,'UserData');
%   pom(Ktory,1:2)=['0' 0];
  pom(Ktory,:)=pom(old_Ktory,:);
  set(hEditO,'UserData',pom);

  pom=get(hEditK,'UserData');
%   pom(Ktory,1:5)=['1024' 0];
  pom(Ktory,:)=pom(old_Ktory,:);
  set(hEditK,'UserData',pom);

  pom=get(hMenu,'UserData');
  pom(Ktory,1:size(pom,2))=ones(1,size(pom,2))*NaN;
  set(hMenu,'UserData',pom);

  pom=get(h_real,'UserData');
%   pom(Ktory,1)=1;
  pom(Ktory,:)=pom(old_Ktory,:);
  set(h_real,'Value', pom(Ktory,1));
  set(h_real,'UserData',pom);

  pom=get(h_det,'UserData');
%   pom(Ktory,1)=0; % \Fixed 2005.11.02
  pom(Ktory,:)=pom(old_Ktory,:);
  set(h_det,'Value', pom(Ktory,1));
  set(h_det,'UserData',pom);

  pom=get(h_rgb(1),'UserData');
  pom(Ktory,1:3)=[0 0 0];
  for ind=1:3,
    set(h_rgb(ind),'Value', pom(Ktory,ind));
  end;
  set(h_rgb(1),'UserData',pom);

  set(hMenu,'Value', Ktory);
  
  % \added 2006.10.09
  pom = get(h_per_type, 'UserData');
  pom(Ktory,:)=pom(old_Ktory,:);
  set(h_per_type,'Value', pom(Ktory,1));
  set(h_per_type,'UserData',pom);
  pom = get(h_Fs, 'UserData');
  pom(Ktory,:)=pom(old_Ktory,:);
  set(h_Fs,'String', pom(Ktory,:));
  set(h_Fs,'UserData',pom);
  % /added 2006.10.09
  
  PERgraf('Choose');
  PERgraf('signal', Ktory);
  return;
  
elseif strcmp(Akcja, 'Reset')
  Ktory=get(hMenu,'Value');  

  pom=get(hMenu,'String');
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
  pom(Ktory,1:5)=['-600' 0];
  set(hEditNoise,'UserData',pom);
  
%  pom=get(hEditN,'UserData');
%  pom(Ktory,1:2)=['1' 0];
%  set(hEditN,'UserData',pom);
  pom=get(hEditO,'UserData');
  pom(Ktory,1:2)=['0' 0];
  set(hEditO,'UserData',pom);

  pom=get(hEditK,'UserData');
  pom(Ktory,1:5)=['1024' 0];
  set(hEditK,'UserData',pom);

  pom=get(hMenu,'UserData');
  for ind=1:3,
    if isfinite(pom(Ktory,ind))
      delete(pom(Ktory,ind));
    end
  end
%   pom(Ktory,:)=[];
  pom(Ktory,1:size(pom,2))=ones(1,size(pom,2))*NaN;
  set(hMenu,'UserData',pom);

  pom=get(h_real,'UserData');
  pom(Ktory,1)=1;
  set(h_real,'Value', pom(Ktory,1));
  set(h_real,'UserData',pom);

  pom=get(h_det,'UserData');
  pom(Ktory,1)=0; % \Fixed 2005.11.02
  set(h_det,'Value', pom(Ktory,1));
  set(h_det,'UserData',pom);

  pom=get(h_rgb(1),'UserData');
  pom(Ktory,1:3)=[0 0 0];
  for ind=1:3,
    set(h_rgb(ind),'Value', pom(Ktory,ind));
  end;
  set(h_rgb(1),'UserData',pom);

  set(hMenu,'Value', Ktory);
  
  % \added 2006.10.09
  pom = get(h_per_type, 'UserData');
  pom(Ktory,1)=1;
  set(h_per_type,'Value', pom(Ktory,1));
  set(h_per_type,'UserData',pom);
  pom = get(h_Fs, 'UserData');
  pom(Ktory,1:4)=['Inf', 0];
  set(h_Fs,'String', pom(Ktory,:));
  set(h_Fs,'UserData',pom);
  % /added 2006.10.09

  PERgraf('Choose');
  PERgraf('signal', Ktory);
  return;
  
elseif strcmp(Akcja, 'delete')
  pom=get(hMenu,'String');
  if size(pom,1)==1,
    PERgraf('Reset');
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
  
  pom=get(h_rgb(1),'UserData');
  pom(Ktory,:)=[];
  set(h_rgb(1),'UserData',pom);
  set(hMenu,'Value', 1);

  % \added 2006.10.09
  pom = get(h_per_type, 'UserData');
  pom(Ktory,:)=[];
  set(h_per_type,'UserData',pom);
  pom = get(h_Fs, 'UserData');
  pom(Ktory,:)=[];
  set(h_Fs,'UserData',pom);
  % /added 2006.10.09
  
  PERgraf('Choose');
  PERgraf('signal', 1); % \Fixed 2005.11.02
  return;
  
% Changing previously entered data set  
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
%   tmp_=['L=' pom(Ktory,:) ';']
%   eval(tmp_, 'L=0;'); % \Fixed 2005.11.02 
  L=eval(pom(Ktory,:), '0'); % \Fixed 2005.11.02 
  
  pom=get(hEditM,'UserData');
  set(hEditM, 'String', pom(Ktory,:));
%   eval(['M=' pom(Ktory,:) ';'], 'M=0;'); % \Fixed 2005.11.02 
  M=eval(pom(Ktory,:), '0'); % \Fixed 2005.11.02 
  
  pom=get(hEditNoise,'UserData');
  set(hEditNoise, 'String', pom(Ktory,:));
  
  pom=get(hEditO,'UserData');
  set(hEditO, 'String', pom(Ktory,:));
%   eval(['O=' pom(Ktory,:) ';'], 'O=0;'); % \Fixed 2005.11.02 
  O=eval(pom(Ktory,:), '0'); % \Fixed 2005.11.02 
  O=round(O/100*M); % \Fixed 2005.11.03 nak³adkowanie podawane w procentach !!! 

%  pom=get(hEditN,'UserData');
%  set(hEditN, 'String', pom(Ktory,:));
  N = fix((L-O)/(M-O)); % \Fixed 2005.11.02 Updates N when user data set is changed
  set(hEditN, 'String', sprintf('%i', N));
  

  pom=get(hEditK,'UserData');
  set(hEditK, 'String', pom(Ktory,:));

  pom=get(h_real,'UserData');
  set(h_real,'Value', pom(Ktory,1));

  pom=get(h_det,'UserData');
  set(h_det,'Value', pom(Ktory,1));

  pom=get(h_rgb(1),'UserData');
  for ind=1:3,
    set(h_rgb(ind),'Value', sign(pom(Ktory,ind)));
  end;

  % \added 2006.10.09
  pom = get(h_per_type, 'UserData');
  per_type = pom(Ktory,1);
  set(h_per_type,'Value', pom(Ktory,1));
  pom = get(h_Fs, 'UserData');
  set(h_Fs,'String', pom(Ktory,:));
  % /added 2006.10.09

  % \Fixed 2005.11.02
  pom=get(hMenu,'UserData');
  psd_mean = pom(:,8); % \Fixed 2005.11.02
  psd_var = pom(:,9); % \Fixed 2005.11.02
  set(hParam,'string',sprintf('mean=%f\r\nvar =%f',psd_mean(Ktory),psd_var(Ktory)));
  pom=get(h_rgb(1),'UserData');
  kolor=pom(Ktory,:);
  set(hParam, 'ForegroundColor', kolor); % \Fixed 2005.11.02
  % /Fixed 2005.11.02
  
%  PERgraf('signal', Ktory);
  UpdatePergraphGUI(per_type, fig);
  return;
  
elseif strcmp(Akcja, 'signal')
  set(fig, 'Pointer', 'watch'); % \Fixed 2005.11.02
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
  hp_per=pom(:,3);
  max_x=pom(:,4);
  min_per=pom(:,5);
  max_per=pom(:,6);
  max_L = pom(:,7); % \Fixed 2005.11.02
  psd_mean = pom(:,8); % \Fixed 2005.11.02
  psd_var = pom(:,9); % \Fixed 2005.11.02

  
  %generate signal
  tekstL=get(hEditL,'UserData');  tekstL=tekstL(Ktory,:);
  ind=find(tekstL==0);
  if ~isempty(ind)
    tekstL=tekstL(1:ind(1)-1);
  end;
%  eval(['L=' tekstL ';'], 'L=1;')
  L = round(getVal(tekstL, 1, hA, 'Error in parameter L')); % \Fixed 2005.11.02
  if L < 1,
    L = 1;
  end
  max_L(Ktory) = L;  % \Fixed 2005.11.02

  A = get(hA.A_edit, 'UserData');
  A = A.A;

  % \added 2006.10.09
  tekst_Fs=get(h_Fs,'UserData');  tekst_Fs=tekst_Fs(Ktory,:);
  ind=find(tekst_Fs==0);
  if ~isempty(ind)
    tekst_Fs=tekst_Fs(1:ind(1)-1);
  end;
  Fs = getVal(tekst_Fs, Inf, hA, ['sampling rate error: >>' tekst_Fs '<<']); 
  % /added 2006.10.09
  
  n=0:L-1;
  tekstX=get(hEditX,'UserData');  tekstX=tekstX(Ktory,:);
  ind=find(tekstX==0);
  if ~isempty(ind)
    tekstX=tekstX(1:ind(1)-1);
  end;
%   eval(['x=' tekstX ';'], ['disp(''signal string error: >>' tekstX '<<''); x=zeros(1,L);']) 
  % \Fixed 2005.11.02
  err_ = 0;
  eval(['x=' tekstX ';'], ['err_=1; x=zeros(1,L);']) 
  if err_==1,
    tmp_ = getVal('NaN', 0, hA, ['signal string error: >>' tekstX '<<']); 
  end
  % /Fixed 2005.11.02
  
  Re=get(h_real,'UserData');  Re=Re(Ktory,1);
  
  tekstNoise=get(hEditNoise,'UserData');  tekstNoise=tekstNoise(Ktory,:);
  ind=find(tekstNoise==0);
  if ~isempty(ind)
    tekstNoise=tekstNoise(1:ind(1)-1);
  end;
%   eval(['Noise=' tekstNoise ';'], 'Noise=-300;') 
  Noise = getVal(tekstNoise, -300, hA, 'Error in parameter Noise'); % \Fixed 2005.11.02

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
    x=x+(N_lin*randn(size(x))+N_lin*j*randn(size(x)))/sqrt(2);
  end;
  max_x(Ktory)=max(abs([real(x); imag(x)]));
   
  %draw signal
  pom=get(h_rgb(1),'UserData');
  kolor=pom(Ktory,:);

  %real part of the signal
  axes(ha(1));
  if isnan(hp_re(Ktory))
    hold on;
    hp_re(Ktory)=plot(n,real(x), 'Color', kolor);
    hold off;
  else
    set(hp_re(Ktory),'Xdata', n, 'Ydata', real(x), 'Color', kolor);
  end
  set(ha(1), 'Xlim', [-0.5, max(max_L)-0.5], 'Ylim', [-1.1 1.1]*max(max_x));  % \Fixed 2005.11.02
  eval('rmappdata(get(ha(1),''Zlabel''),''ZOOMAxesData'')','set(get(ha(1),''ZLabel''),''UserData'',[])');
%  eval('zoom reset', 'set(get(ha(1),''ZLabel''),''UserData'',[]);');    
%  reset(get(ha(1),'ZLabel'));    
  
  %imaginary part of the signal
  axes(ha(2));
  if isnan(hp_im(Ktory))
    hold on;
    hp_im(Ktory)=plot(n,imag(x), 'Color', kolor);
    hold off;
  else
    set(hp_im(Ktory),'Xdata', n, 'Ydata', imag(x), 'Color', kolor);
  end
  set(ha(2), 'Xlim', [-0.5, max(max_L)-0.5], 'Ylim', [-1.1 1.1]*max(max_x));  % \Fixed 2005.11.02
%  set(get(ha(2),'ZLabel'),'UserData',[]);    
%  reset(get(ha(2),'ZLabel'));    
  eval('rmappdata(get(ha(2),''Zlabel''),''ZOOMAxesData'')','set(get(ha(2),''ZLabel''),''UserData'',[])');
%  eval('zoom reset', 'set(get(ha(2),''ZLabel''),''UserData'',[]);');    
  

  set(hMenu,'UserData', [hp_re, hp_im, hp_per, max_x, min_per, max_per, max_L, psd_mean, psd_var]);  % \Fixed 2005.11.02
  
  set(fig, 'Pointer', 'arrow'); % \Fixed 2005.11.02

  %compute and draw periodogram
  PERgraf('per', Ktory)
  return;
  
elseif strcmp(Akcja, 'per_type_change')
  per_type = get(h_per_type, 'Value');
  pom = get(h_per_type, 'UserData');
  pom(Ktory,1)=per_type;
  set(h_per_type,'UserData',pom);

  UpdatePergraphGUI(per_type, fig);
  %compute and draw periodogram
  PERgraf('per', Ktory)
  
elseif strcmp(Akcja, 'per')
  set(fig, 'Pointer', 'watch'); % \Fixed 2005.11.02
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
  hp_per=pom(:,3);
  max_x=pom(:,4);
  min_per=pom(:,5);
  max_per=pom(:,6);
  max_L = pom(:,7);  % \Fixed 2005.11.02
  psd_mean = pom(:,8); % \Fixed 2005.11.02
  psd_var = pom(:,9); % \Fixed 2005.11.02
 
  per_type = get(h_per_type,'UserData');  per_type=per_type(Ktory,1);
  
  %DFT length
  tekstK=get(hEditK,'UserData');  tekstK=tekstK(Ktory,:);
  ind=find(tekstK==0);
  if ~isempty(ind)
    tekstK=tekstK(1:ind(1)-1);
  end;
%   eval(['K=' tekstK ';'], 'K=16;')
  K = getVal(tekstK, 16, hA, 'Error in parameter K'); % \Fixed 2005.11.02

%%%%%%%%%%%%%%%%%%%%%
% \fixed 2006.12.20
%%%%%%%%%%%%%%%%%%%%%
  % signal * BEGIN *
  x=get(hp_re(Ktory), 'Ydata')+j*get(hp_im(Ktory), 'Ydata');
  % signal * END *
  
  if (per_type == 1) | (per_type == 2)
    M = length(x);
    %if (length(x) ~= M)
    %  tmp_ = getVal('NaN', 0, hA, ['signal length is different from periodograph segment length (L != M)']); 
    %end
    %if length(x) < M,
    %  x(M) = 0;
    %else
    %  x = x(1:M);
    %end;
    
    
    if (M>K)
      if supress_msgbox ==1,
        warning('L cannot be larger than K');
      else
        uiwait(msgbox('L cannot be larger than K', 'modal'));
      end
      M=K;
      x = x(1:M);
    end;
    
  else
  
    tekstM=get(hEditM,'UserData');  tekstM=tekstM(Ktory,:);
    ind=find(tekstM==0);
    if ~isempty(ind)
      tekstM=tekstM(1:ind(1)-1);
    end;
  %   eval(['M=' tekstM ';'], 'M=16;') 
    M = round(getVal(tekstM, 16, hA, 'Error in parameter M')); % \Fixed 2005.11.02
    if M < 1,
      M = 1;
    end
    
    if (M>K)
      if supress_msgbox ==1,
        warning('M cannot be larger than K');
      else
        uiwait(msgbox('M cannot be larger than K', 'modal'));
      end
      M=K;
    end;
  end
%%%%%%%%%%%%%%%%%%%%%
% /fixed 2006.12.20
%%%%%%%%%%%%%%%%%%%%%

%  tekstN=get(hEditN,'UserData');  tekstN=tekstN(Ktory,:);
%  ind=find(tekstN==0);
%  if ~isempty(ind)
%    tekstN=tekstN(1:ind(1)-1);
%  end;
%  eval(['N=' tekstN ';'], 'N=16;')

  if per_type == 4,
    tekstO=get(hEditO,'UserData');  tekstO=tekstO(Ktory,:);
    ind=find(tekstO==0);
    if ~isempty(ind)
      tekstO=tekstO(1:ind(1)-1);
    end;
  %   eval(['O=' tekstO ';'], 'O=0;') 
    O = getVal(tekstO, 0, hA, 'Error in parameter O'); % \Fixed 2005.11.02
    O = round(O/100*M); % \Fixed 2005.11.02 nak³adkowanie podawane w procentach !!! 
  else
    O = 0; % no overlapping
  end
  
  if per_type ~= 1,
    tekstW=get(hEditW,'UserData');  tekstW=tekstW(Ktory,:);
    ind=find(tekstW==0);
    if ~isempty(ind)
      tekstW=tekstW(1:ind(1)-1);
    end;

    % \Fixed 2005.11.02
    A = get(hA.A_edit, 'UserData');
    A = A.A;

    err_ = 0;
    eval(['w=' tekstW ';'], 'err_=1; w=ones(1,M);')
    if err_==1,
      tmp_ = getVal('NaN', 0, hA, ['window string error: >>' tekstW '<<']); 
    end
    % /Fixed 2005.11.02
  else
    w = ones(1,M);
  end

  w=w(:);
  if length(w)<M;
    w(M)=0;
  else
    w=w(1:M);
  end;
  
  N = fix((length(x)-O)/(M-O)); 
  set(hEditN, 'String', sprintf('%i', N));

  det_=get(h_det,'UserData');  det_=det_(Ktory,1);
  if det_==0,
    det_='none';
  else
    det_='linear';
  end;
  
  Fs = get(h_Fs, 'UserData'); Fs = Fs(Ktory, :);
  Fs = getVal(Fs, Inf, hA, 'Error in parameter Fs');
  if ~isfinite(Fs),
    Fs = 1;
  end
  %[Per, f]=psd(x, K, Fs, w, O, det_);

  if strcmp(det_, 'none') ~= 1,
    w_tmp = w(:);
    N_tmp = length(x);    % Number of data points
    N_wind = length(w);            % length of window
    % zero-pad x if it has length less than the window length:
    if N_tmp < N_wind            
      x(N_wind)=0;  
      N_tmp = N_wind;
    end
    % Number of windows (k = fix(n/nwind) for noverlap=0):
    k_tmp = fix((N_tmp-O)/(N_wind-O));
    % Zero-pad to allow x to have equal number of elements as k window elements.
    if N_tmp < N_wind*k_tmp,
        x(N_wind*k_tmp)=0;
    end
    xr_tmp=reshape(x,N_wind,k_tmp) % Resizes X column-wise to nwind
    xr_tmp=detrend(xr_tmp);    % DETRENDS each column, aka segment
    x=reshape(xr_tmp,length(x),1) % Makes x back into a column vector
  end
  % Korekta 2018.03.02
  if length(x) < length(w)
      x(length(w)) = 0;
      if supress_msgbox ==1,
        warning('L should not be smaller than M');
      else
        uiwait(msgbox('L should not be smaller than M', 'modal'));
      end
  end
  [Per, f]=pwelch(x, w, O, K, Fs);

  Re=get(h_real,'UserData'); 
  if ~Re(Ktory,1)
    f=fftshift(f);
    f(1:floor(K/2))=f(1:floor(K/2))-Fs;  % \fixed 2006.10.10
    Per=fftshift(Per);
  else
    f=[-f(length(f):-1:1); f];
    Per=[Per(length(Per):-1:1); Per];
  end
  psd_mean(Ktory)=mean(Per);
  psd_var(Ktory)=std(Per).^2;
  set(hParam,'string',sprintf('mean=%f\r\nvar =%f',psd_mean(Ktory),psd_var(Ktory)));
  
  min_per(Ktory)=min([Per(isfinite(Per)); 0]);
  max_per(Ktory)=max([Per(isfinite(Per)); 0.001]);
  if get(h_dB, 'UserData') %dB
    Per=10*log10(Per);
  end
  
  %draw signal
  pom=get(h_rgb(1),'UserData');
  kolor=pom(Ktory,:);
  set(hParam, 'ForegroundColor', kolor); % \Fixed 2005.11.02

  %periodogram
  axes(ha(3));
  if isnan(hp_per(Ktory))
    hold on;
    hp_per(Ktory)=plot(f,Per, 'Color', kolor);
    hold off;
  else
    set(hp_per(Ktory),'Xdata', f, 'Ydata', Per, 'Color', kolor);
  end
  max_Fs = getMaxFs(h_Fs, hA);
  if all(Re)
    set(ha(3), 'Xlim', [0 max_Fs/2]);
  else
    set(ha(3), 'Xlim', [-max_Fs/2 max_Fs/2]);
  end
%  set(get(ha(3),'ZLabel'),'UserData',[]);    
%  reset(get(ha(3),'ZLabel'));    
  eval('rmappdata(get(ha(3),''Zlabel''),''ZOOMAxesData'')','set(get(ha(3),''ZLabel''),''UserData'',[])');
%  eval('zoom reset', 'set(get(ha(3),''ZLabel''),''UserData'',[]);');    
  
  set(hMenu,'UserData', [hp_re, hp_im, hp_per, max_x, min_per, max_per, max_L, psd_mean, psd_var]);  % \Fixed 2005.11.02

  PERgraf('per_ylim');
  zoom on;
  set(fig, 'Pointer', 'arrow'); % \Fixed 2005.11.02
  return;
  
elseif strcmp(Akcja, 'Color')
  pom=get(h_rgb(1),'UserData');
  for ind=1:3,
    pom(Ktory,ind)=get(h_rgb(ind),'Value')*0.75;
  end;
  set(h_rgb(1),'UserData', pom);
    
  pom2=get(hMenu,'UserData');
  hp_re=pom2(:,1);
  hp_im=pom2(:,2);
  hp_per=pom2(:,3);

  %PERgraf('redraw', Ktory); %w zasadzie tutaj tylko zmiana koloru bez przeliczania
  set([hp_re(Ktory), hp_im(Ktory), hp_per(Ktory)], 'Color', pom(Ktory,:));
  set(hParam, 'ForegroundColor', pom(Ktory,:)); % \Fixed 2005.11.02
  return;
  
elseif strcmp(Akcja, 'lock')
  lock = get(h_lock, 'Value');
  if (lock == 1)
    pom = [get(ha(3), 'Xlim'), get(ha(3), 'Ylim')];
    set(h_lock, 'UserData', pom);
  end
  PERgraf('per_ylim');
  
elseif strcmp(Akcja, 'dB')
  pom=get(h_dB, 'UserData');
    if pom~=get(h_dB, 'Value')
    pom=get(hMenu,'UserData');
    hp_per=pom(:,3);
    for ind=1:length(hp_per)
      sygn=get(hp_per(ind),'Ydata');
      if get(h_dB, 'Value') %dB
        sygn=10*log10(sygn);
      else %lin
        sygn=10.^(sygn/10);
      end;
      set(hp_per(ind),'Ydata', sygn);
    end;
    set(h_dB, 'UserData', get(h_dB, 'Value'));
    
    lock = get(h_lock, 'Value');
    if (lock == 1)
      pom= get(h_lock, 'UserData');
      if get(h_dB, 'Value') %dB
        pom(3:4) = 10*log10(pom(3:4));
        if ~isfinite(pom(3)),
          pom(3) = -300;
        end
      else %lin
        pom(3:4) = 10.^(pom(3:4)/10);
      end;
      set(h_lock, 'UserData', pom);
    end
  end
  PERgraf('per_ylim');
  return
  
elseif strcmp(Akcja, 'per_ylim')
  lock = get(h_lock, 'Value');
  if (lock == 1)
    pom= get(h_lock, 'UserData');
    if length(pom) ~= 4,
      set(h_lock, 'Value', 0);
      PERgraf('lock');
      return;
    end
    set(ha(3), 'Xlim', pom(1:2));
    set(ha(3), 'Ylim', pom(3:4));
  else
    pom=get(hMenu,'UserData');
    min_per=pom(:,5);
    max_per=pom(:,6);
    if get(h_dB, 'UserData') %dB
      ylim_=[min(min_per) max(max_per)*1.1];
      if isfinite(ylim_(2))
        ylim_(2)=2.^ceil(log2(ylim_(2)));
      else
        ylim_(2)=0;
      end
      warning off;
      ylim_=10*log10(ylim_);
      warning on;
      tekst=get(hEdit_dY,'String');
      eval(['dY=' tekst ';'], 'dY=120;'); 
      if dY<=0, dY=10; end;
      if isfinite(ylim_(1))
        ylim_(1)=2.^floor(log2(ylim_(1)));
      else
        ylim_(1)=ylim_(2)-dY;
      end
      if ylim_(2)-ylim_(1)>dY,
        ylim_(1)=ylim_(2)-dY;
      end
      set(ha(3), 'Ylim', ylim_);
    else
      ylim_=2.^ceil(log2(max(max_per)));
      set(ha(3), 'Ylim', [0 ylim_]);
    end
    
    Re=get(h_real,'UserData'); 
    max_Fs = getMaxFs(h_Fs, hA);
    if all(Re)
      set(ha(3), 'Xlim', [0 max_Fs/2]);
    else
      set(ha(3), 'Xlim', [-max_Fs/2 max_Fs/2]);
    end
  end
  pause(0);
  % eval('zoom reset', 'set(get(ha(3),''ZLabel''),''UserData'',[]);');
  eval('rmappdata(get(ha(3),''Zlabel''),''ZOOMAxesData'')','set(get(ha(3),''ZLabel''),''UserData'',[])');
  % zoom on;
  return
  
elseif strcmp(Akcja, 'actionA')
  redraw_ = 0;
  A = get(hA.A_edit, 'UserData');
  if strcmp(param, 'A_edit')
    %   - "A_edit"
    A_ = eval(get(hA.A_edit, 'String'), 'NaN');
    if ~isfinite(A_),
      set(hA.A_edit, 'String', A.A);
    else
      A.A = A_;
      if A.A < A.min;
        A.A = A.min;
        set(hA.A_edit, 'String', A.A);
      end
      if A.A > A.max;
        A.A = A.max;
        set(hA.A_edit, 'String', A.A);
      end
      set(hA.A_slider, 'Value', A.A);
      redraw_ = 1;
    end

  elseif strcmp(param, 'dA_edit')
    %   - "dA_edit"
    dA_ = eval(get(hA.dA_edit, 'String'), 'NaN');
    if ~isfinite(dA_),
      set(hA.dA_edit, 'String', A.dA);
    else
      A.dA = dA_;
    end
    
  elseif strcmp(param, 'A_slider')
    %   - "A_slider"
    A_ = get(hA.A_slider, 'Value');
    znak = sign(A_ - A.A);
    if A.is_plus == 1,
      A_ = A.A + znak*A.dA;        
    else % razy
      if znak >= 0,
        A_ = A.A*A.dA;        
      else
        A_ = A.A/A.dA;        
      end
    end

    A.A = A_;
    if A.A < A.min;
      A.A = A.min;
    end
    if A.A > A.max;
      A.A = A.max;
    end
    set(hA.A_edit, 'String', A.A);
    set(hA.A_slider, 'Value', A.A);
    redraw_ = 1;
    
  elseif strcmp(param, 'minA_edit')
    %   - "minA_edit"
    minA_ = eval(get(hA.minA_edit, 'String'), 'NaN');
    if ~isfinite(minA_),
      set(hA.minA_edit, 'String', A.min);
    else
      if minA_ > A.max;
        minA_ = A.min;
        set(hA.minA_edit, 'String', minA_);
      end
      
      A.min = minA_;
      set(hA.A_slider, 'min', minA_);
      
      if A.A < A.min;
        A.A = A.min;
        set(hA.A_edit, 'String', A.A);
        set(hA.A_slider, 'Value', A.A);
      end
    end
    
  elseif strcmp(param, 'maxA_edit')
    %   - "maxA_edit"
    maxA_ = eval(get(hA.maxA_edit, 'String'), 'NaN');
    if ~isfinite(maxA_),
      set(hA.maxA_edit, 'String', A.max);
    else
      if maxA_ < A.min;
        maxA_ = A.max;
        set(hA.maxA_edit, 'String', maxA_);
      end
      
      A.max = maxA_;
      set(hA.A_slider, 'max', maxA_);
      
      if A.A > A.max;
        A.A = A.max;
        set(hA.A_edit, 'String', A.A);
        set(hA.A_slider, 'Value', A.A);
      end
    end
    
  elseif strcmp(param, 'A_razy_radio')
    %   - "A_razy_radio"
%     val =  get(A_plus_radio, 'Value');
%     if val > 0.0,
%       A.is_plus = 1;
%     else
      A.is_plus = 0;
%     end
%     if A.is_plus == 1,
%       set(hA.A_plus_radio, 'Value', 1.0);
%       set(hA.A_razy_radio, 'Value', 0.0);
%     else
      set(hA.A_plus_radio, 'Value', 0.0);
      set(hA.A_razy_radio, 'Value', 1.0);
%     end
    
  elseif strcmp(param, 'A_plus_radio')
    %   - "A_plus_radio"
%     val =  get(A_razy_radio, 'Value');
%     if val > 0.0,
%       A.is_plus = 0;
%     else
      A.is_plus = 1;
%     end
%     if A.is_plus == 1,
      set(hA.A_plus_radio, 'Value', 1.0);
      set(hA.A_razy_radio, 'Value', 0.0);
%     else
%       set(hA.A_plus_radio, 'Value', 0.0);
%       set(hA.A_razy_radio, 'Value', 1.0);
%     end
    
  end
  
  set(hA.A_edit, 'UserData', A);
  supress_msgbox = 1;
  if redraw_ == 1,
    Ktory=get(hMenu,'Value');  
    pom=get(hMenu,'UserData');
    ile=size(pom,1);
    for ind=1:ile,
      if ind ~= Ktory,
        eval(sprintf('pergraf(''signal'', %i);', ind), 'set(fig, ''Pointer'', ''arrow'');');
      end
    end
    eval(sprintf('pergraf(''signal'', %i);', Ktory), 'set(fig, ''Pointer'', ''arrow'');');
  end
  supress_msgbox = 0;
  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = getVal(tekst, def_val, hA, msg, Fs)

global supress_msgbox

if nargin == 4,
  Fs = Inf;
end

A = get(hA.A_edit, 'UserData');
A = A.A;

x = eval(tekst, 'NaN');
if isnan(x),
  x=def_val;
  if supress_msgbox ==1,
    warning(msg);
  else
    uiwait(msgbox(msg, 'modal'));
  end
end

function UpdatePergraphGUI(per_type, fig)

hEditM=findobj(fig,'tag', 'M_edit');
hLabelM=findobj(fig,'tag', 'M_label');
hEditN=findobj(fig,'tag', 'N_edit');
hLabelN=findobj(fig,'tag', 'N_label');
hEditO=findobj(fig,'tag', 'O_edit');
hLabelO=findobj(fig,'tag', 'O_label');
hUnitO=findobj(fig,'tag', 'O_unit');
hEditW=findobj(fig,'tag', 'w_edit');
hLabelW=findobj(fig,'tag', 'w_label');

switch (per_type)
  case 1, % raw periodograph
    set(hEditM, 'Visible', 'off');
    set(hLabelM, 'Visible', 'off');
    set(hEditN, 'Visible', 'off');
    set(hLabelN, 'Visible', 'off');
    set(hEditO, 'Visible', 'off');
    set(hLabelO, 'Visible', 'off');
    set(hUnitO, 'Visible', 'off');
    set(hEditW, 'Visible', 'off');
    set(hLabelW, 'Visible', 'off');
  case 2, % modified periodograph
    set(hEditM, 'Visible', 'off');
    set(hLabelM, 'Visible', 'off');
    set(hEditN, 'Visible', 'off');
    set(hLabelN, 'Visible', 'off');
    set(hEditO, 'Visible', 'off');
    set(hLabelO, 'Visible', 'off');
    set(hUnitO, 'Visible', 'off');
    set(hEditW, 'Visible', 'on');
    set(hLabelW, 'Visible', 'on');
  case 3, % Bartlett periodograph
    set(hEditM, 'Visible', 'on');
    set(hLabelM, 'Visible', 'on');
    set(hEditN, 'Visible', 'on');
    set(hLabelN, 'Visible', 'on');
    set(hEditO, 'Visible', 'off');
    set(hLabelO, 'Visible', 'off');
    set(hUnitO, 'Visible', 'off');
    set(hEditW, 'Visible', 'on');
    set(hLabelW, 'Visible', 'on');
  case 4, % Welch periodograph
    set(hEditM, 'Visible', 'on');
    set(hLabelM, 'Visible', 'on');
    set(hEditN, 'Visible', 'on');
    set(hLabelN, 'Visible', 'on');
    set(hEditO, 'Visible', 'on');
    set(hLabelO, 'Visible', 'on');
    set(hUnitO, 'Visible', 'on');
    set(hEditW, 'Visible', 'on');
    set(hLabelW, 'Visible', 'on');
  otherwise,
    disp('Unsupported periodograph type');
end;

function maxFs = getMaxFs(h_Fs, hA)

maxFs = 0;
pom = get(h_Fs, 'UserData');
for ind=1:size(pom, 1);
  tekst_Fs = pom(ind, :);
  ind=find(tekst_Fs==0);
  if ~isempty(ind)
    tekst_Fs = tekst_Fs(1:ind(1)-1);
  end;
  Fs = getVal(tekst_Fs, Inf, hA, ['sampling rate error: >>' tekst_Fs '<<']); 
  if ~isfinite(Fs),
    Fs = 1;
  end;
  maxFs = max([maxFs, Fs]);
end;
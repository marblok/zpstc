% Wielostopniowa zmiana szybkoœci próbkowania
function MultiSRC_gui_v2
 
fig = findobj('tag', 'MultiSRC_gui_v2b');
if ~isempty(fig)
   figure(fig);
   return;
end
clc
set(0, 'DefaultUicontrolFontName', 'Times New Roman CE')

global fold
global fnew
 
global L
global M
 
global ile_stopni 
 
global nr_sygnalu
 
 
function result = is_calkowita(liczba)
 
if liczba == fix(liczba)
 result = true;
 else
 result = false;
 end
end
 
function keyPress(hObject, eventdata)
    switch eventdata.Key
    	case 'f11'
            guihelp;
    end
end


function Help(hObject, eventdata)
    guihelp;
end
 
 function Odswiez(hObject, eventdata)
 
 fold = str2double(get(edtFold, 'String'));
 fnew = str2double(get(edtFnew, 'String'));
 
 if isnan(fold) || isnan(fnew) || ~is_calkowita(fold) || ~is_calkowita(fnew)
 msgbox('Nieprawid³owa wartoœæ Fold lub Fnew.', 'B³¹d');
 return;
 end
 
 L = fnew / gcd(fold, fnew);
 M = fold / gcd(fold, fnew); 

 Fmax = min(0.5/L,0.5/M)*(fold*L);
 set(txtFmax, 'String', sprintf('F_max = %.1f [Hz]', Fmax));
 
 set(txtL, 'String', sprintf('L = %d', L));
 set(txtM, 'String', sprintf('M = %d', M));
 
 set(edtL1, 'String', '');
 set(edtL2, 'String', '');
 set(edtL3, 'String', '');
 set(edtL4, 'String', '');
 set(edtM1, 'String', '');
 set(edtM2, 'String', '');
 set(edtM3, 'String', '');
 set(edtM4, 'String', '');
 set(edtFP1, 'String', '');
 set(edtFP2, 'String', '');
 set(edtFP3, 'String', '');
 set(edtFP4, 'String', '');
 set(edtFS1, 'String', '');
 set(edtFS2, 'String', '');
 set(edtFS3, 'String', '');
 set(edtFS4, 'String', '');
%  for i = 2:7
%     if ishandle(i)
%         close(i);
%     end
%  end
 end


function WyborSygnalu(hObject, eventdata)
    nr_sygnalu = get(hObject, 'value');
end

function WyliczDomyslneFsFp(hObject, eventdata)

Fp_factor = str2double(get(edtFp_factor, 'String'));

L1 = str2double(get(edtL1, 'String'));
M1 = str2double(get(edtM1, 'String'));
L2 = str2double(get(edtL2, 'String')); 
M2 = str2double(get(edtM2, 'String'));
L3 = str2double(get(edtL3, 'String')); 
M3 = str2double(get(edtM3, 'String'));
L4 = str2double(get(edtL4, 'String')); 
M4 = str2double(get(edtM4, 'String'));

if (ile_stopni == 1)
    Li = L1;
    Mi = M1;    
    
    if isnan(L1) || isnan(M1) || ~is_calkowita(L1) || ~is_calkowita(M1) 
        msgbox('Nieprawid³owa wartoœæ w Li lub Mi.', 'B³¹d');
        return;
    end
    
    if prod(Li) ~= L || prod(Mi) ~= M
        msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
        return;
    end
    
    % Wartoœci referencyjne implementacji jednostopniowej
%     Fp = min(0.4/L,0.4/M)*(fold*L1);
    Fs = min(0.5/L,0.5/M)*(fold*L); % filtr dla implementacji jednostopniowej pracuje na czêstotliwoœci Fold*L
    Fp = (1-Fp_factor)*Fs;
    Fmax = Fs; % maksymalna czêstotliwoœæ sk³adowych sygna³u na wyjœciu implementacji jednostopniowej
    
    % Wartoœci dla pierwszego stopnia
% 	Fp1 = min(0.4/L1,0.4/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fs1 = min(0.5/L1,0.5/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fp1 = (1-Fp_factor)*Fs1;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs1 > Fmax),
        Fs1 = Fs1 + (Fs1-Fmax);
    end
    if (Fp1 >= Fp),
      Fp1 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp1 < Fp');
    end
    if (Fs1 < Fmax)
        disp('Fs1 < Fmax');
    end
    fp1 = Fp1/(fold*L1); fs1 = Fs1/(fold*L1);
    Fmax1 = min(Fmax, Fs1);
    if (Fmax1 < Fmax)
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax1), 'ForegroundColor', 'red', 'FontWeight', 'Bold');
    else
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax1), 'ForegroundColor', 'black', 'FontWeight', 'normal');
    end
    
    set(edtFP1, 'String', fp1);
    set(edtFS1, 'String', fs1);
end
 
if (ile_stopni == 2)
    Li = [ L1, L2 ];
    Mi = [ M1, M2 ];

    if isnan(L1) || isnan(M1) || ~is_calkowita(L1) || ~is_calkowita(M1) ||...
       isnan(L2) || isnan(M2) || ~is_calkowita(L2) || ~is_calkowita(M2)            
        msgbox('Nieprawid³owa wartoœæ w Li lub Mi.', 'B³¹d');
        return;
    end
    
    if prod(Li) ~= L || prod(Mi) ~= M
        msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
        return;
    end
    
    % Wartoœci referencyjne implementacji jednostopniowej
%     Fp = min(0.4/L,0.4/M)*(fold*L1);
    Fs = min(0.5/L,0.5/M)*(fold*L); % filtr dla implementacji jednostopniowej pracuje na czêstotliwoœci Fold*L
%     Fp = 0.8*Fs;
    Fp = (1-Fp_factor)*Fs;
    Fmax = Fs; % maksymalna czêstotliwoœæ sk³adowych sygna³u na wyjœciu implementacji jednostopniowej
    
    % Wartoœci dla pierwszego stopnia
% 	Fp1 = min(0.4/L1,0.4/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fs1 = min(0.5/L1,0.5/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fp1 = (1-Fp_factor)*Fs1;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs1 >= Fmax),
        Fs1 = Fs1 + (Fs1-Fmax);
    end
    if (Fp1 >= Fp),
      Fp1 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp1 < Fp');
    end
    if (Fs1 < Fmax)
        disp('Fs1 < Fmax');
    end
    fp1 = Fp1/(fold*L1); fs1 = Fs1/(fold*L1);
    Fmax1 = min(Fmax, Fs1);
    if (Fmax1 < Fmax)
        disp('Fmax1 < Fmax');
    end
    
    
    % Wartoœci dla drugiego stopnia
% 	Fp2 = min(0.4/L2,0.4/M2)*(fold*L1/M1*L2); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci fold*L1/M1*L2
    Fs2 = min(0.5/L2,0.5/M2)*(fold*L1/M1*L2); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1/M1*L2
    Fp2 = (1-Fp_factor)*Fs2;
    % alias jest dopuszczalny do czêstotliwoœci Fmax
    if (Fs2 >= Fmax),
        Fs2 = Fs2 + (Fs2-Fmax);
    end
    if (Fp2 >= Fp),
      Fp2 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp2 < Fp');
    end
    if (Fs2 < Fmax)
        disp('Fs2 < Fmax');
    end
    fp2 = Fp2/(fold*L1/M1*L2); fs2 = Fs2/(fold*L1/M1*L2);
    Fmax2 = min(Fmax1, Fs2);
    if (Fmax2 < Fmax)
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax2), 'ForegroundColor', 'red', 'FontWeight', 'Bold');
    else
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax2), 'ForegroundColor', 'black', 'FontWeight', 'normal');
    end
    
% fp1 = (min(0.4/L,0.4/M)*fnew)/(fold*L1);
% fs1 = min(0.5/L1,0.5/M1);
% fp2 = (min(0.4/L,0.4/M)*fnew)/(fold*L1*L2/M1);
% fs2 = min(0.5/L2,0.5/M2);
    set(edtFP1, 'String', fp1);
    set(edtFS1, 'String', fs1);
    set(edtFP2, 'String', fp2);
    set(edtFS2, 'String', fs2);
end
 
if (ile_stopni == 3)
    Li = [ L1, L2, L3 ];
    Mi = [ M1, M2, M3 ];

    if isnan(L1) || isnan(M1) || ~is_calkowita(L1) || ~is_calkowita(M1) ||...
       isnan(L2) || isnan(M2) || ~is_calkowita(L2) || ~is_calkowita(M2) ||...            
       isnan(L3) || isnan(M3) || ~is_calkowita(L3) || ~is_calkowita(M3)
        msgbox('Nieprawid³owa wartoœæ w Li lub Mi.', 'B³¹d');
        return;
    end
    
    if prod(Li) ~= L || prod(Mi) ~= M
        msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
        return;
    end
    
    % Wartoœci referencyjne implementacji jednostopniowej
%     Fp = min(0.4/L,0.4/M)*(fold*L1);
    Fs = min(0.5/L,0.5/M)*(fold*L); % filtr dla implementacji jednostopniowej pracuje na czêstotliwoœci Fold*L
    Fp = (1-Fp_factor)*Fs;
    Fmax = Fs; % maksymalna czêstotliwoœæ sk³adowych sygna³u na wyjœciu implementacji jednostopniowej
    
%    fp1 = (min(0.4/L,0.4/M)*fnew)/(fold*L1);
% fs1 = min(0.5/L1,0.5/M1);
% 	Fp1 = min(0.4/L1,0.4/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fs1 = min(0.5/L1,0.5/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fp1 = (1-Fp_factor)*Fs1;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs1 >= Fmax),
        Fs1 = Fs1 + (Fs1-Fmax);
    end
    if (Fp1 >= Fp),
      Fp1 = Fp; % wybieramy szersze pasmo przejœciowe
    else
     disp('Fp1 < Fp');
    end
    if (Fs1 < Fmax)
        disp('Fs1 < Fmax');
    end
    fp1 = Fp1/(fold*L1); fs1 = Fs1/(fold*L1);
    Fmax1 = min(Fmax, Fs1);
    if (Fmax1 < Fmax)
        disp('Fmax1 < Fmax');
    end
    
% fp2 = (min(0.4/L,0.4/M)*fnew)/(fold*L1*L2/M1);
% fs2 = min(0.5/L2,0.5/M2);
    % Wartoœci dla drugiego stopnia
% 	Fp2 = min(0.4/L2,0.4/M2)*(fold*L1/M1*L2); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci fold*L1/M1*L2
    Fs2 = min(0.5/L2,0.5/M2)*(fold*L1/M1*L2); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1/M1*L2
    Fp2 = (1-Fp_factor)*Fs2;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs2 >= Fmax),
        Fs2 = Fs2 + (Fs2-Fmax);
    end
    if (Fp2 >= Fp),
      Fp2 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp2 < Fp');
    end
    if (Fs2 < Fmax)
        disp('Fs2 < Fmax');
    end
    fp2 = Fp2/(fold*L1/M1*L2); fs2 = Fs2/(fold*L1/M1*L2);
    Fmax2 = min(Fmax1, Fs2);
    if (Fmax2 < Fmax)
        disp('Fmax2 < Fmax');
    end
    
% fp3 = (min(0.4/L,0.4/M)*fnew)/(fold*L1*L2*L3/(M1*M2));
% fs3 = min(0.5/L3,0.5/M3);
    % Wartoœci dla TRZECIEGO stopnia
% 	Fp3 = min(0.4/L3,0.4/M3)*(fold*L1/M1*L2/M2*L3); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci fold*L1/M1*L2
    Fs3 = min(0.5/L3,0.5/M3)*(fold*L1/M1*L2/M2*L3); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1/M1*L2
    Fp3 = (1-Fp_factor)*Fs3;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs3 >= Fmax),
        Fs3 = Fs3 + (Fs3-Fmax);
    end
    if (Fp3 >= Fp),
      Fp3 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp3 < Fp');
    end
    if (Fs3 < Fmax)
        disp('Fs3 < Fmax');
    end
    fp3 = Fp3/(fold*L1/M1*L2/M2*L3); fs3 = Fs3/(fold*L1/M1*L2/M2*L3);
    Fmax3 = min(Fmax2, Fs3);
    if (Fmax3 < Fmax)
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax3), 'ForegroundColor', 'red', 'FontWeight', 'Bold');
    else
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax3), 'ForegroundColor', 'black', 'FontWeight', 'normal');
    end


    set(edtFP1, 'String', fp1);
    set(edtFS1, 'String', fs1);
    set(edtFP2, 'String', fp2);
    set(edtFS2, 'String', fs2);
    set(edtFP3, 'String', fp3);
    set(edtFS3, 'String', fs3);
end
 
if (ile_stopni == 4)
    Li = [ L1, L2, L3, L4 ];
    Mi = [ M1, M2, M3, M4 ];

    if isnan(L1) || isnan(M1) || ~is_calkowita(L1) || ~is_calkowita(M1) ||...
       isnan(L2) || isnan(M2) || ~is_calkowita(L2) || ~is_calkowita(M2) ||...            
       isnan(L3) || isnan(M3) || ~is_calkowita(L3) || ~is_calkowita(M3) ||...
       isnan(L4) || isnan(M4) || ~is_calkowita(L4) || ~is_calkowita(M4)
        msgbox('Nieprawid³owa wartoœæ w Li lub Mi.', 'B³¹d');
        return;
    end
    
    if prod(Li) ~= L || prod(Mi) ~= M
        msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
        return;
    end
    
    % Wartoœci referencyjne implementacji jednostopniowej
%     Fp = min(0.4/L,0.4/M)*(fold*L1);
    Fs = min(0.5/L,0.5/M)*(fold*L); % filtr dla implementacji jednostopniowej pracuje na czêstotliwoœci Fold*L
    Fp = (1-Fp_factor)*Fs;
    Fmax = Fs; % maksymalna czêstotliwoœæ sk³adowych sygna³u na wyjœciu implementacji jednostopniowej
    
% fp1 = (min(0.4/L,0.4/M)*fnew)/(fold*L1);
% fs1 = min(0.5/L1,0.5/M1);
% 	Fp1 = min(0.4/L1,0.4/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fs1 = min(0.5/L1,0.5/M1)*(fold*L1); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1
    Fp1 = (1-Fp_factor)*Fs1;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs1 >= Fmax),
        Fs1 = Fs1 + (Fs1-Fmax);
    end
    if (Fp1 >= Fp),
      Fp1 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp1 < Fp');
    end
    if (Fs1 < Fmax)
        disp('Fs1 < Fmax');
    end
    fp1 = Fp1/(fold*L1); fs1 = Fs1/(fold*L1);
    Fmax1 = min(Fmax, Fs1);
    if (Fmax1 < Fmax)
        disp('Fmax1 < Fmax');
    end
    
% fp2 = (min(0.4/L,0.4/M)*fnew)/(fold*L1*L2/M1);
% fs2 = min(0.5/L2,0.5/M2);
    % Wartoœci dla drugiego stopnia
% 	Fp2 = min(0.4/L2,0.4/M2)*(fold*L1/M1*L2); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci fold*L1/M1*L2
    Fs2 = min(0.5/L2,0.5/M2)*(fold*L1/M1*L2); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1/M1*L2
    Fp2 = (1-Fp_factor)*Fs2;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs2 >= Fmax),
        Fs2 = Fs2 + (Fs2-Fmax);
    end
    if (Fp2 >= Fp),
      Fp2 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp2 < Fp');
    end
    if (Fs2 < Fmax)
        disp('Fs2 < Fmax');
    end
    fp2 = Fp2/(fold*L1/M1*L2); fs2 = Fs2/(fold*L1/M1*L2);
    Fmax2 = min(Fmax1, Fs2);
    if (Fmax2 < Fmax)
        disp('Fmax2 < Fmax');
    end
    
% fp3 = (min(0.4/L,0.4/M)*fnew)/(fold*L1*L2*L3/(M1*M2));
% fs3 = min(0.5/L3,0.5/M3);
    % Wartoœci dla TRZECIEGO stopnia
% 	Fp3 = min(0.4/L3,0.4/M3)*(fold*L1/M1*L2/M2*L3); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci fold*L1/M1*L2
    Fs3 = min(0.5/L3,0.5/M3)*(fold*L1/M1*L2/M2*L3); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1/M1*L2
    Fp3 = (1-Fp_factor)*Fs3;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs3 >= Fmax),
        Fs3 = Fs3 + (Fs3-Fmax);
    end
    if (Fp3 >= Fp),
      Fp3 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp3 < Fp');
    end
    if (Fs3 < Fmax)
        disp('Fs3 < Fmax');
    end
    fp3 = Fp3/(fold*L1/M1*L2/M2*L3); fs3 = Fs3/(fold*L1/M1*L2/M2*L3);
    Fmax3 = min(Fmax2, Fs3);
    if (Fmax3 < Fmax)
        disp('Fmax3 < Fmax');
    end
    
% fp4 = (min(0.4/L,0.4/M)*fnew)/(fold*L1*L2*L3*L4/(M1*M2*M3))
% fs4 = (min(0.5/L,0.5/M)*fnew)/(fold*L1*L2*L3*L4/(M1*M2*M3))
    % Wartoœci dla czwartego stopnia
% 	Fp4 = min(0.4/L4,0.4/M4)*(fold*L1/M1*L2/M2*L3/M3*L4); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci fold*L1/M1*L2
    Fs4 = min(0.5/L4,0.5/M4)*(fold*L1/M1*L2/M2*L3/M3*L4); % filtr dla pierwszego stopnia pracuje na czêstotliwoœci Fold*L1/M1*L2
    Fp4 = (1-Fp_factor)*Fs4;
    % alias jest dopuszczalny do szybkoœci Fmax
    if (Fs4 >= Fmax),
        Fs4 = Fs4 + (Fs4-Fmax);
    end
    if (Fp4 >= Fp),
      Fp4 = Fp; % wybieramy szersze pasmo przejœciowe
    else
      disp('Fp4 < Fp');
    end
    if (Fs4 < Fmax)
        disp('Fs4 < Fmax');
    end
    fp4 = Fp4/(fold*L1/M1*L2/M2*L3/M3*L4); fs4 = Fs4/(fold*L1/M1*L2/M2*L3/M3*L4);
    Fmax4 = min(Fmax3, Fs4);
    if (Fmax4 < Fmax)
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax4), 'ForegroundColor', 'red', 'FontWeight', 'Bold');
    else
      set(txtFmax2, 'String', sprintf('F_max = %.1f [Hz]', Fmax4), 'ForegroundColor', 'black', 'FontWeight', 'normal');
    end
    
    set(edtFP1, 'String', fp1);
    set(edtFS1, 'String', fs1);
    set(edtFP2, 'String', fp2);
    set(edtFS2, 'String', fs2);
    set(edtFP3, 'String', fp3);
    set(edtFS3, 'String', fs3);
    set(edtFP4, 'String', fp4);
    set(edtFS4, 'String', fs4);
end

end

 
function WyborStopni(hObject, eventdata)
 
set(edtL1, 'Enable', 'off');
set(edtM1, 'Enable', 'off');
set(edtFP1, 'Enable', 'off');
set(edtFS1, 'Enable', 'off');
 
set(edtL2, 'Enable', 'off');
set(edtM2, 'Enable', 'off');
set(edtFP2, 'Enable', 'off');
set(edtFS2, 'Enable', 'off');
 
set(edtL3, 'Enable', 'off');
set(edtM3, 'Enable', 'off');
set(edtFP3, 'Enable', 'off');
set(edtFS3, 'Enable', 'off');
 
set(edtL4, 'Enable', 'off');
set(edtM4, 'Enable', 'off');
set(edtFP4, 'Enable', 'off');
set(edtFS4, 'Enable', 'off'); 
 
set(edtL1, 'String', '');
set(edtM1, 'String', '');
set(edtFP1, 'String', '');
set(edtFS1, 'String', '');
 
set(edtL2, 'String', '');
set(edtM2, 'String', '');
set(edtFP2, 'String', '');
set(edtFS2, 'String', '');
 
set(edtL3, 'String', '');
set(edtM3, 'String', '');
set(edtFP3, 'String', '');
set(edtFS3, 'String', '');
 
set(edtL4, 'String', '');
set(edtM4, 'String', '');
set(edtFP4, 'String', '');
set(edtFS4, 'String', '');
 
ile_stopni = get(hObject, 'Value');
 
if (ile_stopni > 0)
set(edtL1, 'Enable', 'on');
set(edtM1, 'Enable', 'on');
set(edtFP1, 'Enable', 'on');
set(edtFS1, 'Enable', 'on');
end
 
if (ile_stopni > 1)
set(edtL2, 'Enable', 'on');
set(edtM2, 'Enable', 'on');
set(edtFP2, 'Enable', 'on');
set(edtFS2, 'Enable', 'on');
 end
 
 if (ile_stopni > 2)
 set(edtL3, 'Enable', 'on');
 set(edtM3, 'Enable', 'on');
 set(edtFP3, 'Enable', 'on');
 set(edtFS3, 'Enable', 'on');
 end
 
 if (ile_stopni > 3)
 set(edtL4, 'Enable', 'on');
 set(edtM4, 'Enable', 'on');
 set(edtFP4, 'Enable', 'on');
 set(edtFS4, 'Enable', 'on');
 end
 
 end
 

 
 function Nnum = ObliczNnum(Ni, Li, Mi)
 stopien = length(Ni);
 
 if (stopien >= 1), FLPFi(1) = fold * Li(1); end
 if (stopien >= 2), FLPFi(2) = fold * Li(1) * Li(2) / Mi(1); end
 if (stopien >= 3), FLPFi(3) = fold * Li(1) * Li(2) * Li(3) / (Mi(1) * Mi(2)); end
 if (stopien >= 4), FLPFi(4) = fold * Li(1) * Li(2) * Li(3) * Li(4) / (Mi(1) * Mi(2) * Mi(3)); end
 
 Nnum = sum(Ni .* FLPFi) / fnew;
 Nnum = ceil(Nnum);
 end
 

 
 function Start(hObject, eventdata)
 
 set(txtN1, 'String', '');
 set(txtN2, 'String', '');
 set(txtN3, 'String', '');
 set(txtN4, 'String', '');
 set(txtNtot, 'String', '');
 set(txtNnum, 'String', '');
 drawnow;
 
%  for i = 2:7
%     if ishandle(i)
%         close(i);
%     end
%  end
 
 if (ile_stopni == 1)
 L1 = str2double(get(edtL1, 'String')); 
 M1 = str2double(get(edtM1, 'String'));
 
 FP1 = str2double(get(edtFP1, 'String'));
 FS1 = str2double(get(edtFS1, 'String'));
 
 Li = L1;
 Mi = M1;
 
 FPi = FP1;
 FSi = FS1;
 
 if isnan(L1) || isnan(M1) || isnan(FP1) || isnan(FS1) || ~is_calkowita(L1) || ~is_calkowita(M1) 
 msgbox('Nieprawid³owa wartoœæ w Li, Mi, FPi lub FSi.', 'B³¹d');
 return;
 end
 
 if FP1 >= FS1
 msgbox('FSi musi byc wiêksze od FPi', 'B³¹d');
 return;
 end
 
 if prod(Li) ~= L || prod(Mi) ~= M
 msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
 return;
 end
 
 N = id_1(nr_sygnalu, fold, fnew, Li, Mi, FPi, FSi);
 set(txtN1, 'String', sprintf('N1 = %d', N(1)));
 end
 
 if (ile_stopni == 2)
 L1 = str2double(get(edtL1, 'String')); 
 M1 = str2double(get(edtM1, 'String'));
 M2 = str2double(get(edtM2, 'String'));
 L2 = str2double(get(edtL2, 'String')); 
 
 FP1 = str2double(get(edtFP1, 'String'));
 FS1 = str2double(get(edtFS1, 'String'));
 FP2 = str2double(get(edtFP2, 'String'));
 FS2 = str2double(get(edtFS2, 'String'));
 
 Li = [ L1, L2 ];
 Mi = [ M1, M2 ];
 
 FPi = [ FP1, FP2 ];
 FSi = [ FS1, FS2 ];
 
 if isnan(L1) || isnan(M1) || isnan(FP1) || isnan(FS1) || ~is_calkowita(L1) || ~is_calkowita(M1) ||...
 isnan(L2) || isnan(M2) || isnan(FP2) || isnan(FS2) || ~is_calkowita(L2) || ~is_calkowita(M2)
 msgbox('Nieprawid³owa wartoœæ w Li, Mi, FPi lub FSi.', 'B³¹d');
 return;
 end
 
 if FP1 >= FS1 || FP2 >= FS2
 msgbox('FSi musi byc wiêksze od FPi', 'B³¹d');
 return;
 end
 
 if prod(Li) ~= L || prod(Mi) ~= M
 msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
 return;
 end
 
 N = id_2(nr_sygnalu, fold, fnew, Li, Mi, FPi, FSi);
 
 set(txtN1, 'String', sprintf('N1 = %d', N(1)));
 set(txtN2, 'String', sprintf('N2 = %d', N(2)));
 end
 
 if (ile_stopni == 3)
 L1 = str2double(get(edtL1, 'String')); 
 M1 = str2double(get(edtM1, 'String'));
 L2 = str2double(get(edtL2, 'String')); 
 M2 = str2double(get(edtM2, 'String'));
 L3 = str2double(get(edtL3, 'String')); 
 M3 = str2double(get(edtM3, 'String'));
 
 FP1 = str2double(get(edtFP1, 'String'));
 FS1 = str2double(get(edtFS1, 'String'));
 FP2 = str2double(get(edtFP2, 'String'));
 FS2 = str2double(get(edtFS2, 'String'));
 FP3 = str2double(get(edtFP3, 'String'));
 FS3 = str2double(get(edtFS3, 'String'));
 
 Li = [ L1, L2, L3 ];
 Mi = [ M1, M2, M3 ];
 
 FPi = [ FP1, FP2, FP3 ];
 FSi = [ FS1, FS2, FS3 ];
 
 if isnan(L1) || isnan(M1) || isnan(FP1) || isnan(FS1) || ~is_calkowita(L1) || ~is_calkowita(M1) ||...
 isnan(L2) || isnan(M2) || isnan(FP2) || isnan(FS2) || ~is_calkowita(L2) || ~is_calkowita(M2) ||...
 isnan(L3) || isnan(M3) || isnan(FP3) || isnan(FS3) || ~is_calkowita(L3) || ~is_calkowita(M3)
 msgbox('Nieprawid³owa wartoœæ w Li, Mi, FPi lub FSi.', 'B³¹d');
 return;
 end
 
 if FP1 >= FS1 || FP2 >= FS2 || FP3 >= FS3
 msgbox('FSi musi byc wiêksze od FPi', 'B³¹d');
 return;
 end
 
 if prod(Li) ~= L || prod(Mi) ~= M
 msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
 return;
 end
 
 N = id_3(nr_sygnalu, fold, fnew, Li, Mi, FPi, FSi);
 
 set(txtN1, 'String', sprintf('N1 = %d', N(1)));
 set(txtN2, 'String', sprintf('N2 = %d', N(2)));
 set(txtN3, 'String', sprintf('N3 = %d', N(3)));
 end 
 
 if (ile_stopni == 4)
 L1 = str2double(get(edtL1, 'String')); 
 M1 = str2double(get(edtM1, 'String'));
 L2 = str2double(get(edtL2, 'String')); 
 M2 = str2double(get(edtM2, 'String'));
 L3 = str2double(get(edtL3, 'String')); 
 M3 = str2double(get(edtM3, 'String'));
 L4 = str2double(get(edtL4, 'String')); 
 M4 = str2double(get(edtM4, 'String'));
 
 FP1 = str2double(get(edtFP1, 'String'));
 FS1 = str2double(get(edtFS1, 'String'));
 FP2 = str2double(get(edtFP2, 'String'));
 FS2 = str2double(get(edtFS2, 'String'));
 FP3 = str2double(get(edtFP3, 'String'));
 FS3 = str2double(get(edtFS3, 'String'));
 FP4 = str2double(get(edtFP4, 'String'));
 FS4 = str2double(get(edtFS4, 'String'));
 
 Li = [ L1, L2, L3, L4 ];
 Mi = [ M1, M2, M3, M4 ];
 
 FPi = [ FP1, FP2, FP3, FP4 ];
 FSi = [ FS1, FS2, FS3, FS4 ];
 
 if isnan(L1) || isnan(M1) || isnan(FP1) || isnan(FS1) || ~is_calkowita(L1) || ~is_calkowita(M1) ||...
 isnan(L2) || isnan(M2) || isnan(FP2) || isnan(FS2) || ~is_calkowita(L2) || ~is_calkowita(M2) ||...
 isnan(L3) || isnan(M3) || isnan(FP3) || isnan(FS3) || ~is_calkowita(L3) || ~is_calkowita(M3) ||... 
 isnan(L4) || isnan(M4) || isnan(FP4) || isnan(FS4) || ~is_calkowita(L4) || ~is_calkowita(M4)
 msgbox('Nieprawid³owa wartoœæ w Li, Mi, FPi lub FSi.', 'B³¹d');
 return;
 end
 
 if FP1 >= FS1 || FP2 >= FS2 || FP3 >= FS3 || FP4 >= FS4
 msgbox('FSi musi byc wiêksze od FPi', 'B³¹d');
 return;
 end
 
 if prod(Li) ~= L || prod(Mi) ~= M
 msgbox('Iloczyn Li (Mi) musi byæ równy L (M).', 'B³¹d');
 return;
 end
 
 N = id_4(nr_sygnalu, fold, fnew, Li, Mi, FPi, FSi);
 
 set(txtN1, 'String', sprintf('N1 = %d', N(1)));
 set(txtN2, 'String', sprintf('N2 = %d', N(2)));
 set(txtN3, 'String', sprintf('N3 = %d', N(3)));
 set(txtN4, 'String', sprintf('N4 = %d', N(4))); 
 end
  
 Ntot = sum(N);
 Nnum = ObliczNnum(N, Li, Mi);
 
 set(txtNtot, 'String', sprintf('Ntot = %d', Ntot));
 set(txtNnum, 'String', sprintf('Nnum = %d', Nnum));
 
 end
 
function onCloseRequest(x, y)
  yn = questdlg('Zamkn¹æ aplikacjê?',...
      'Figure Menu',...
      'Yes','No','No');
  switch yn
      case 'Yes'
        main_fig = findobj('tag', 'MultiSRC_gui_v2b');
        if ~isempty(main_fig)
          delete(main_fig);
        end
      case 'No'
        return
  end
end

function onDelete(~,~)
    main_fig = findobj('tag', 'MultiSRC_gui_v2b');
    if ~isempty(main_fig)
      set(main_fig,'DeleteFcn',[]);
      delete(main_fig);
    end
        
    for ind = 2:7
      subfig = findobj('tag', sprintf("MultiSRC_gui_v2b_fig%i", ind));
      if ~isempty(subfig)
        % set(subfig,'DeleteFcn',[]);
        delete(subfig);
      end
    end
end
 
screen_resolution = get(0, 'ScreenSize');
screen_width = screen_resolution(3);
screen_height = screen_resolution(4);
 
window_width = 427;
window_height = 525+40;
 
window_position = [(screen_width/2 - window_width/2),...
 (screen_height/2 - window_height/2),...
 window_width,...
 window_height];
 
main = figure('Name','Konwersja szybkoœci próbkowania (ver. 2b)',...
 'NumberTitle','off',...
 'HandleVisibility', 'on', ...
 'IntegerHandle', 'off', ...
 'Resize','off',...
 'Visible','on',...
 'MenuBar', 'none',...
 'DoubleBuffer', 'on',...
 'units', 'pixels',...
 'color', [.9 .9 .9],...
 'KeyPressFcn', @keyPress,...
 'DeleteFcn', @onDelete,...
 'CloseRequestFcn', @onCloseRequest,...
 'Position', window_position,...
 'tag', 'MultiSRC_gui_v2b');
 

 
panel1 = uipanel('Parent', main, 'Units', 'pixels', 'Visible', 'on', 'Position', [5 320+40 420 200]);
 
txtL = uicontrol('Parent', panel1, 'Style', 'text', 'String', 'L = ', 'Position', [10 35 150 30], 'FontSize', 12, 'FontWeight', 'Bold');
txtM = uicontrol('Parent', panel1, 'Style', 'text', 'String', 'M = ', 'Position', [210 35 150 30], 'FontSize', 12, 'FontWeight', 'Bold');
txtFmax = uicontrol('Parent', panel1, 'Style', 'text', 'String', 'F_max = ', 'Position', [70 15 170 20], 'FontSize', 12, 'HorizontalAlignment', 'left');

uicontrol('Parent', panel1, 'Style', 'pushbutton', 'String', 'Odœwie¿', 'FontName', 'Times New Roman CE', 'Position', [120 70 150 30], 'Callback', @Odswiez);
 
uicontrol('Parent', panel1, 'Style', 'text', 'String', 'F_in =  [Sa/S]', 'Position', [10 120 150 15], 'HorizontalAlignment', 'left');
edtFold = uicontrol('Parent', panel1, 'Style', 'edit', 'String', '', 'Position', [85 115 65 25], 'BackgroundColor', 'white', 'HorizontalAlignment', 'left');
 
uicontrol('Parent', panel1, 'Style', 'text', 'String', 'F_out = [Sa/S]', 'Position', [210 120 150 15], 'HorizontalAlignment', 'left');
edtFnew = uicontrol('Parent', panel1, 'Style', 'edit', 'String', '', 'Position', [295 115 65 25], 'BackgroundColor', 'white', 'HorizontalAlignment', 'left');
 
uicontrol('Parent', panel1, 'Style', 'text', 'String', 'SYGNA£ TESTOWY', 'FontName', 'Times New Roman CE', 'Position', [80 160 100 15], 'HorizontalAlignment', 'left');
uicontrol('Parent', panel1, 'Style', 'popupmenu', 'String', {'Sygna³ 1', 'Sygna³ 2'}, 'FontName', 'Times New Roman CE', 'Position', [190 160 100 20], 'Callback', @WyborSygnalu, 'BackgroundColor', 'white');
 
uicontrol('Parent', panel1, 'Style', 'pushbutton', 'String', 'POMOC  [F11]', 'Position', [330 155 80 30], 'Callback', @Help);
 
panel2 = uipanel('Parent', main, 'Units', 'pixels', 'Visible', 'on', 'Position', [5 155 420 160+40]);
 
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'L4 = ', 'Position', [10 15+40 25 15], 'HorizontalAlignment', 'left');
edtL4 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [40 10+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'M4 = ', 'Position', [110 15+40 25 15], 'HorizontalAlignment', 'left');
edtM4 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [140 10+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
uicontrol('Parent', panel2, 'Style', 'text', 'String', '<fp,fs> LPF4   <                  ,                   >', 'Position', [210 15+40 200 15], 'HorizontalAlignment', 'left');
edtFP4 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [290 10+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
edtFS4 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [350 10+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
 
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'L3 = ', 'Position', [10 40+40 25 15], 'HorizontalAlignment', 'left');
edtL3 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [40 35+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'M3 = ', 'Position', [110 40+40 25 15], 'HorizontalAlignment', 'left');
edtM3 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [140 35+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
uicontrol('Parent', panel2, 'Style', 'text', 'String', '<fp,fs> LPF3   <                  ,                   >', 'Position', [210 40+40 200 15], 'HorizontalAlignment', 'left');
edtFP3 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [290 35+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
edtFS3 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [350 35+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
 
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'L2 = ', 'Position', [10 65+40 25 15], 'HorizontalAlignment', 'left');
edtL2 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [40 60+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'M2 = ', 'Position', [110 65+40 25 15], 'HorizontalAlignment', 'left');
edtM2 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [140 60+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
uicontrol('Parent', panel2, 'Style', 'text', 'String', '<fp,fs> LPF2   <                  ,                   >', 'Position', [210 65+40 200 15], 'HorizontalAlignment', 'left');
edtFP2 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [290 60+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
edtFS2 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [350 60+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'off');
 
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'L1 = ', 'Position', [10 90+40 25 15], 'HorizontalAlignment', 'left');
edtL1 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [40 85+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'on');
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'M1 = ', 'Position', [110 90+40 25 15], 'HorizontalAlignment', 'left');
edtM1 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [140 85+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'on');
uicontrol('Parent', panel2, 'Style', 'text', 'String', '<fp,fs> LPF1   <                  ,                   >', 'Position', [210 90+40 200 15], 'HorizontalAlignment', 'left');
edtFP1 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [290 85+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'on');
edtFS1 = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '', 'Position', [350 85+40 50 25], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white', 'Enable', 'on');
 
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'LICZBA STOPNI', 'Position', [15 130+40 100 15], 'HorizontalAlignment', 'left');
uicontrol('Parent', panel2, 'Style', 'popupmenu', 'String', {'1', '2', '3', '4'}, 'Position', [100 130+40 50 20], 'Callback', @WyborStopni, 'BackgroundColor', 'white');
uicontrol('Parent', panel2, 'Style', 'text', 'String', 'alfa = ', 'Position', [160 130+40 40 15], 'HorizontalAlignment', 'right');
edtFp_factor = uicontrol('Parent', panel2, 'Style', 'edit', 'String', '0.2', 'Position', [205 125+40 50 25], 'BackgroundColor', 'white', 'HorizontalAlignment', 'left');
uicontrol('Parent', panel2, 'Style', 'pushbutton', 'String', 'Wylicz fp i fs', 'Position', [290 122+40 110 30], 'Callback', @WyliczDomyslneFsFp);

txtFmax2 = uicontrol('Parent', panel2, 'Style', 'text', 'String', '', 'Position', [70 15 170 20], 'FontSize', 12, 'HorizontalAlignment', 'left');

 
panel3 = uipanel('Parent', main, 'Units', 'pixels', 'Visible', 'on', 'Position', [5 5 420 145]);
 
txtN1 = uicontrol('Parent', panel3, 'Style', 'text', 'String', '', 'Position', [10 50 100 30], 'FontSize', 12, 'FontWeight', 'Bold');
txtN2 = uicontrol('Parent', panel3, 'Style', 'text', 'String', '', 'Position', [110 50 100 30], 'FontSize', 12, 'FontWeight', 'Bold');
txtN3 = uicontrol('Parent', panel3, 'Style', 'text', 'String', '', 'Position', [210 50 100 30], 'FontSize', 12, 'FontWeight', 'Bold');
txtN4 = uicontrol('Parent', panel3, 'Style', 'text', 'String', '', 'Position', [310 50 100 30], 'FontSize', 12, 'FontWeight', 'Bold');
 
txtNtot = uicontrol('Parent', panel3, 'Style', 'text', 'String', '', 'Position', [ 35 10 150 30], 'FontSize', 12, 'FontWeight', 'Bold');
txtNnum = uicontrol('Parent', panel3, 'Style', 'text', 'String', '', 'Position', [235 10 150 30], 'FontSize', 12, 'FontWeight', 'Bold');
 
uicontrol('Parent', panel3, 'Style', 'pushbutton', 'String', 'START', 'Position', [50 100 100 30], 'Callback', @Start);
uicontrol('Parent', panel3, 'Style', 'text', 'String', 'uruchom zmianê szybkoœci próbkowania', 'FontName', 'Times New Roman CE', 'Position', [160 100 250 20], 'HorizontalAlignment', 'left');
 

 
nr_sygnalu = 1;
ile_stopni = 1;

end


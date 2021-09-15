function test_symbol_mapper(offset)
if nargin == 0,
  offset = 0;
end

[symb, Fp] = fileread('../symbols.flt', 'cplx');
% % [symb, Fp] = fileread('../symbols.flt'); % dla ASK
% [symb, Fp] = fileread('../ex6_task2_symb_a.flt', 'cplx');

% [y{3}, Fp] = fileread('../ex6_task2_symb_c.flt', 'cplx');

h = fopen('../Ex6_task1.cpp');
bits_in = fread(h, 80, 'ubit1', 'ieee-be').'; % MSB first
bits_in2 = fread(h, 80, 'ubit1', 'ieee-le').'; % LSB first
fclose(h);

% !!! kolejno�� wczytywania bit�w jest odmienna od standardowej!!!
h = fopen('../decoded.txt');
% h = fopen('../ex6_task2_bina.txt');
bits_out = fread(h, 80, 'ubit1', 'ieee-be').'; % MSB first // here the bits offset is corrected in the *.cpp
fclose(h);
h = fopen('../decoded_raw.txt');
% h = fopen('../ex6_task2_bina.txt');
bits_raw = fread(h, 80, 'ubit1', 'ieee-be').'; % MSB first // here the bits offset is corrected in the *.cpp
fclose(h);
h = fopen('../bits_out.txt');
bits_out2 = fread(h, 80, 'char').';
fclose(h);
bits_out2 = double(bits_out2-'0') % For SymbolDemapper this gives 2 additional bits at the begining 
% h = fopen('../bits3.flt');
% bits_out3 = fread(h, 80, 'float').';
% fclose(h);

figure(1)
subplot(1,1,1)
plot(symb, 'o')

bits_out_mod = bits_out;
bits_out_mod(1:2:end) = bits_out(2:2:end);
bits_out_mod(2:2:end) = bits_out(1:2:end);
% bits_out wyprzedza o 2bity (bits_per_symbol) strumienie bits_raw/bits_out2 
figure(2)
subplot(4,1,1)
stem(bits_in, 'b')
subplot(4,1,2)
stem(bits_in2, 'rx')
subplot(4,1,3)
stem(bits_out, 'r')
subplot(4,1,4)
stem(bits_out_mod, 'g')
% subplot(3,1,3)
% stem(bits_raw, 'mx')
% hold on
% stem(bits_out2, 'b+')
% hold off
% % subplot(4,1,4)
% % stem(bits_out3, 'g')

offset1 = 1
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% okre�l przyporz�dkowanie bit�w do symboli w mapperze
bits_ = bits_out2;
% bits_ = bits_out;
if offset1 > 0
  y_ref = symb((offset1+1):end);
else
  y_ref = symb;
end
if offset1 < 0
  if rem(-offset1, 2) == 0,
    bits_ = bits_in((-offset1+1):end);
  else
    bits_ = bits_in((-offset1+1):end-1);
  end
else
  bits_ = bits_in;
end
%bin_ref_new = 2*bits_in-1;
bin_ref_new = bits_;

% \TODO wydoby� konstelacj� wej�ciow� i por�wna� j� z konstelacj� wyj�ciow�
symb_ind = 2*bin_ref_new(1:2:end) + bin_ref_new(2:2:end);
figure(23)
subplot(1,2,1)
kolor = 'rmbg';
scale = [0.9, 0.92, 0.94, 0.96];
for ind = 0:3
  u{ind+1} = scale(ind+1) * y_ref(find(symb_ind == ind));
  plot(real(u{ind+1}), imag(u{ind+1}), [kolor(ind+1),'o'])
  hold on
end
hold off
legend({'00', '10', '01', '11'})
axis equal
title('mapper constellation')

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% okre�l przyporz�dkowanie bit�w do symboli w demapperze
bits_ = bits_out2;
% bits_ = bits_out;
if offset > 0
  y_ref = symb((offset+1):end);
else
  y_ref = symb;
end
if offset < 0
  if rem(-offset, 2) == 0,
    bits_ = bits_((-offset+1):end);
  else
    bits_ = bits_((-offset+1):end-1);
  end
else
  bits_ = bits_;
end
%bin_ref_new = 2*bits_in-1;
bin_ref_new = bits_;

% pierwszy - starszy bit / drugi - m�odszy bit
% % symb_ind = 2*(bin_ref_new(1:2:end)+1)/2 + (bin_ref_new(2:2:end)+1)/2;
symb_ind = 2*bin_ref_new(1:2:end) + bin_ref_new(2:2:end);
figure(23)
subplot(1,2,2)
kolor = 'rmbg';
scale = [0.9, 0.92, 0.94, 0.96];
for ind = 0:3
  u{ind+1} = scale(ind+1) * y_ref(find(symb_ind == ind));
  plot(real(u{ind+1}), imag(u{ind+1}), [kolor(ind+1),'o'])
  hold on
end
hold off
legend({'00', '10', '01', '11'})
axis equal
title('demapper constellation')


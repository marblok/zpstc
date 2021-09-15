function test_task2_constelation(offset)
% test konstelacji dla wszystkich podkana��w
% * test_task2_constelation(2) - dla demodulacji ex6_task1.flt  [mode = 'A']
% * test_task2_constelation(3) - dla demodulacji ex6_task1b.flt [mode = 'B']
if nargin == 0,
  offset = 0;
end

mode = 'A' % demodulacja ex6_task1.flt
% mode = 'B' % demodulacja ex6_task1b.flt

% Surowe symbole z poszczeg�lnych wyodr�bnionych podkanałów Ex6_task_zad2.cpp
[y{1}, Fp] = fileread(['../ex6_task2_symb_', mode, '_a.flt'], 'cplx');
[y{2}, Fp] = fileread(['../ex6_task2_symb_', mode, '_b.flt'], 'cplx');
[y{3}, Fp] = fileread(['../ex6_task2_symb_', mode, '_c.flt'], 'cplx');

% symbole wyj�ciowe mappera symboli dla pierwszego podkanałów Ex6_taskzad1.cpp
%x{1} = fileread('../ex6_task1_ch1.flt', 'cplx');

% symbole wyj�ciowe mappera symboli dla pierwszego podkanałów Ex6_taskzad1b.cpp
x{1} = fileread('../ex6_task1a_new.flt', 'cplx');

y{1} = y{1}((9+1):end);

bin_ref = fileread('../ex6_task1_bin_a.flt', 'cplx');

h = fopen('../ex6_task2_bina.txt');
bits_out = fread(h, inf, 'ubit1', 'ieee-be').';
fclose(h);
h = fopen('../ex6_task2_bina.flt');
bits_out2 = fread(h, inf, 'float').';
fclose(h);

N = length(y{1});
M = length(x{1});

figure(1)
subplot(2,3,1)
plot(y{1}, 'b.')
axis equal
subplot(2,3,2)
plot(y{2}, 'b.')
axis equal
subplot(2,3,3)
plot(y{3}, 'b.')
axis equal

subplot(2,3,4)
plot(x{1}, 'b.')
axis equal

figure(2)
subplot(2,1,1)
plot(real(y{1}), 'b')
hold on
plot(real(x{1}), 'r')
hold off
subplot(2,1,2)
plot(imag(y{1}), 'b')
hold on
plot(imag(x{1}), 'r')
hold off


N1 = min([N,M]);
figure(3)
subplot(2,3,1)
plot(y{1}.*conj(x{1}(1:N1)), 'r.')
axis equal
subplot(2,1,2)
[c, l] = xcorr(y{1}, x{1}(1:N1));
plot(l, real(c), 'b');
hold on
plot(l, imag(c), 'r');
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
symb = y{1};
bits_ = bits_out;
if offset > 0
  y_ref = symb((offset+1):end);
else
  y_ref = symb;
end
if offset < 0
  if rem(-offset, 2) == 0,
    bits_in = bits_((-offset+1):end);
  else
    bits_in = bits_((-offset+1):end-1);
  end
else
  bits_in = bits_;
end
%bin_ref_new = 2*bits_in-1;
bin_ref_new = bits_in;

if length(bin_ref_new) > 2*length(y_ref),
  bin_ref_new = bin_ref_new(1:2*length(y_ref));
end
% pierwszy - starszy bit / drugi - m�odszy bit
% % symb_ind = 2*(bin_ref_new(1:2:end)+1)/2 + (bin_ref_new(2:2:end)+1)/2;
symb_ind = 2*bin_ref_new(1:2:end) + bin_ref_new(2:2:end);
figure(3)
subplot(2,3,3)
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

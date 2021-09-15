function test_task1b_new

x{1} = fileread('../ex6_task1a.flt', 'cplx');
x_new{1} = fileread('../ex6_task1a_new.flt', 'cplx');

[bin_ref, Fp] = fileread('../ex6_task1_bin_a.flt', 'cplx');
[bin_ref_new, Fp] = fileread('../ex6_task1_bin_a_new.flt'); % real
B = length(bin_ref_new);
if rem(B,2) == 1,
  disp('rem(B,2) == 1');
  bin_ref_new = bin_ref_new(1:end-1);
end

figure(1)
subplot(2,2,1)
plot(x{1}, 'bo')
axis equal
subplot(2,2,2)
plot(x_new{1}, 'ro')
axis equal
title('Wej�ciowe strumienie symboli dla kana�u nr 1')

figure(2)
subplot(2,1,1)
plot(real(x{1}), 'b')
hold on
plot(real(x_new{1}), 'r')
hold off
set(gca, 'Xlim', [0,100], 'Ylim', [-1.1, 1.1])
subplot(2,1,2)
plot(imag(x{1}), 'b')
hold on
plot(imag(x_new{1}), 'r')
hold off
set(gca, 'Xlim', [0,100], 'Ylim', [-1.1, 1.1])

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% okre�l przyporz�dkowanie bit�w do symboli modulatorze
y_ref = x{1};
% real - starszy bit / imag - m�odszy bit
symb_ind = (2*real(bin_ref)+1)/2 + (imag(bin_ref)+1)/2;
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
legend({'00', '01', '10', '11'})
axis equal

figure(23)
subplot(1,2,2)
% (2+8+1+1)*2 - wyprzedzenie bin data w bitach
% 
% bin_ref_new = bin_ref_new(3:end);
% +++++++++++++++++++++++++++++++++++++++++++ 
% Bez skip w kodzie 2 bity z przedo za du�o
prefix_size = 2;
bin_ref_new = [-ones(prefix_size,1); bin_ref_new(1:end-prefix_size)];
% % Przy skip w kodzie powinno by� (2+8+1+1)*2-2 bit�w za ma�o - nie ten kod
% prefix_size = (2+8+1)*2;
% bin_ref_new = bin_ref_new(prefix_size+1:end);
y_ref = x_new{1};
% % bin_ref_new = bin_ref_new(2:end-1);
% bin_ref_new = bin_ref_new(3:end);

% pierwszy - starszy bit / drugi - m�odszy bit
symb_ind = 2*(bin_ref_new(1:2:end)+1)/2 + (bin_ref_new(2:2:end)+1)/2;
kolor = 'rmbg';
scale = [0.9, 0.92, 0.94, 0.96];
hp = [];
for ind = 0:3
  u{ind+1} = scale(ind+1) * y_ref(find(symb_ind == ind));
  hp(ind+1)=plot(real(u{ind+1}), imag(u{ind+1}), [kolor(ind+1),'o']);
  hold on
end
hold off
legend(hp, {'00', '01', '10', '11'})
axis equal


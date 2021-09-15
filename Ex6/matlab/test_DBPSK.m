function test_DBPSK

x = readaudiofile('../ex6_task2_text_in.flt');
x_en = readaudiofile('../ex6_task2_text_en.flt', 'cplx');
x_de = readaudiofile('../ex6_task2_text_de.flt');

figure(1)
subplot(3,1,1)
stem(x(1:20))
subplot(3,1,2)
stem(real(x_en(1:20)))
hold on
stem(imag(x_en(1:20)), 'r')
hold off

subplot(3,1,3)
stem(2*(x_de(1:20)-0.5))

return
y(1) = 1;
for ind = 1:length(x_en),
  y(ind+1) = y(1)*x_en(ind);
end

x_en2 = [1; x_en];
x_en2d = [x_en; 1];

y = x_en2 .* conj(x_en2d);

hold on
stem(real(y(1:20))*0.9, 'r')
hold off

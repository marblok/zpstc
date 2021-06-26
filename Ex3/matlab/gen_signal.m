function x = gen_signal(nr_sygnalu, fold, L, M)

N = 5*fold;
while (N < 1024) | (N*L/M < 1024),
    N = N*2;
end

n = 0:N-1;
switch nr_sygnalu
    case 1,
        x = sin(2*pi*9/fold*n)+1/4*sin((2*pi*27/fold*n)+pi/4)+1/16*sin((2*pi*36/fold*n)+pi/6)+sin(2*pi*49/fold*n);
        x = x + randn(size(x))/100000; 
    case 2,
		s = 6;
		x = 0;
        while s < fold
            a = randn;
			x = x + a*sin(2*pi*s/fold*n);
            s = s + 6;
        end
         x = x + randn(size(x))/100000; 
    case 3,
        % x = 
    case 4,
        % x = 
end

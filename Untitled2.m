y = zeros(3, 10)
x = 0.1:0.1:1
y(1, :) = exp(x)
y(2, :) = sin(x)
y(3, :) = 5
semilogy(x, y)
grid on
legend('config 1', 'config 2', 'config 3')

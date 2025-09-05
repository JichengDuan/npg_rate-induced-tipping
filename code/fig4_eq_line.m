
r = 2.0;
K = 3.0;
m = 0.02;
b =1.5;              
% dX/dt - X ͼ
X = linspace(0, 3, 1000);
dX = r .* X .* (1 - X ./ K) .* (b - X) - m .* X;
figure;
plot(X, dX, 'k', 'LineWidth', 2); hold on;


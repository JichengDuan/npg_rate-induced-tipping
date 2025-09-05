syms x
r = 2.0;
K = 3.0;
b = 1.5;
m = 0.02;

f = @(X) r * (1 - X / K) .* (b - X) - m;

a = r / K;
b_coef = -r * (1 + b / K);
c = r * b - m;


Xroots = roots([a, b_coef, c]);


X_eq = [0; Xroots];

X_eq = X_eq(imag(X_eq) == 0); % 实数解
X_eq = real(X_eq);
X_eq = X_eq(X_eq >= 0); 

% 找稳态
fprintf('找到 %d 个稳态：\n', length(X_eq));
for i = 1:length(X_eq)
    Xs = X_eq(i);
    % 计算导数稳定性判断
    dfdX = @(X) r*(1 - X/K)*( -1 ) + r*(-1/K)*(b - X) + r*(1 - X/K)*(-1) - m;
    dF = dfdX(Xs);
    if dF < 0
        stability = 'stable';
    elseif dF > 0
        stability = 'unstalbe';
    else
        stability = 'critical';
    end
    
    fprintf(' X = %.4f ：%s （dF/dX = %.4f）\n', Xs, stability, dF);
end

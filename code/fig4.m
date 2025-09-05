clc;
clear all;
% paramaters
r = 2.0;
K = 3.0;
m = 0.02;
b_start = 1.5;
b_end = 2.0;
Tr = 1.0;
Tend = 10;
v = (b_end - b_start) / Tr;

dt = 0.001;
tspan = 0:dt:Tend;
Nt = length(tspan);

%X0 
X0 = linspace(0, 3, 200);
NX = length(X0);
X_all = zeros(Nt, NX);
%%Rate
for j = 1:NX
    X = X0(j);
    for i = 1:Nt
        t = tspan(i);
        if t <= Tr
            b = b_start + v * t;
        else
            b = b_end;
        end
        dX = r * X * (1 - X / K) * (b - X) - m * X;
        X = X + dt * dX;
        X_all(i, j) = X;
    end
end

X2_track = nan(1, Nt);
X3_track = nan(1, Nt);
for i = 1:Nt
    t = tspan(i);
    if t <= Tr
        b = b_start + v * t;
    else
        b = b_end;
    end
    Delta = (K + b)^2 / 4 - (b * K + m * K / r); %root judge
    if Delta >= 0
        sqrtD = sqrt(Delta);
        X2 = (K + b) / 2 - sqrtD;
        X3 = (K + b) / 2 + sqrtD;
        X2_track(i) = X2;
        X3_track(i) = X3;
    end
end

% 冷冻
b_frozen = b_start;
Delta0 = (K + b_frozen)^2 / 4 - (b_frozen * K + m * K / r);
saddle0 = (K + b_frozen) / 2 - sqrt(Delta0);
% 判定：根据 X0 与初始鞍点比较
colors = repmat("b", 1, NX);
colors(X0 < saddle0) = "r";



figure;
hold on
for j = 1:NX
    if colors(j) == "r"
        plot(tspan, X_all(:, j), 'r');
    else
        plot(tspan, X_all(:, j), 'b');
    end
end
plot(tspan, X2_track, 'm', 'LineWidth', 2);  % 鞍点
plot(tspan, X3_track, 'c', 'LineWidth', 2);  % 高稳态X3
xlabel('time t')
ylabel('X')


clc; clear;

% ---------------------- 参数设置 -------------------------
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon = 1.0;
c_start = 0.4;
c_end = 0.1;
T_r = 8.0;       % 参数漂移时间
T_end = 200;     % 总积分时间

% ---------------------- 初值选择 -------------------------
ICs = [
    0, 4;     % 吸引子右上方
    0, 2.5;     % nullcline 附近
    1.7, 1.5      % 另一个 basin 区域
    1.5, 0
];

% ---------------------- 仿真设置 -------------------------
tspan = [0 T_end];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

% ---------------------- 定义系统 -------------------------
sys = @(t, x) [
    - (x(1)-s1)*(x(1)-s2)*(x(1)-s3) + c_t(t,c_start,c_end,T_r)*x(2);
    - epsilon * (x(2)-s1)*(x(2)-s2)*(x(2)-s3)
];

% ---------------------- 绘图准备 -------------------------
figure; hold on;

% 起始和结束 nullcline 曲线
x1 = linspace(0, 4, 100);
f_start = ((x1 - s1).*(x1 - s2).*(x1 - s3)) / c_start;
f_end   = ((x1 - s1).*(x1 - s2).*(x1 - s3)) / c_end;

plot(x1, f_start, 'b', 'LineWidth', 1.5);         % 起始 nullcline
plot(x1, f_end, 'g', 'LineWidth', 1.5);           % 结束 nullcline

% X2-nullclines
% for x2_nc = [s1, s2, s3]
%     yline(x2_nc, 'k--', 'LineWidth', 1.0);
% end

% ---------------------- 轨道模拟 -------------------------
colors = lines(size(ICs,1));
for i = 1:size(ICs,1)
    x0 = ICs(i,:);
    [t, x] = ode45(sys, tspan, x0, options);
    plot(x(:,1), x(:,2), 'k', 'LineWidth', 1.5);
    plot(x0(1), x0(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
end

% ---------------------- 图形修饰 -------------------------
xlabel('$X_1$', 'Interpreter','latex');
ylabel('$X_2$', 'Interpreter','latex');

axis([0 4 0 4]);
set(gca,'FontSize',12);
box on
function c = c_t(t, c0, c1, Tr)
    if t <= Tr
        c = c0 + (c1 - c0) * (t / Tr);
    else
        c = c1;
    end
end

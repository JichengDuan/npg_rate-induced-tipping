clc; clear;

% ---------------------- �������� -------------------------
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon = 1.0;
c_start = 0.4;
c_end = 0.1;
T_r = 8.0;       % ����Ư��ʱ��
T_end = 200;     % �ܻ���ʱ��

% ---------------------- ��ֵѡ�� -------------------------
ICs = [
    0, 4;     % ���������Ϸ�
    0, 2.5;     % nullcline ����
    1.7, 1.5      % ��һ�� basin ����
    1.5, 0
];

% ---------------------- �������� -------------------------
tspan = [0 T_end];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

% ---------------------- ����ϵͳ -------------------------
sys = @(t, x) [
    - (x(1)-s1)*(x(1)-s2)*(x(1)-s3) + c_t(t,c_start,c_end,T_r)*x(2);
    - epsilon * (x(2)-s1)*(x(2)-s2)*(x(2)-s3)
];

% ---------------------- ��ͼ׼�� -------------------------
figure; hold on;

% ��ʼ�ͽ��� nullcline ����
x1 = linspace(0, 4, 100);
f_start = ((x1 - s1).*(x1 - s2).*(x1 - s3)) / c_start;
f_end   = ((x1 - s1).*(x1 - s2).*(x1 - s3)) / c_end;

plot(x1, f_start, 'b', 'LineWidth', 1.5);         % ��ʼ nullcline
plot(x1, f_end, 'g', 'LineWidth', 1.5);           % ���� nullcline

% X2-nullclines
% for x2_nc = [s1, s2, s3]
%     yline(x2_nc, 'k--', 'LineWidth', 1.0);
% end

% ---------------------- ���ģ�� -------------------------
colors = lines(size(ICs,1));
for i = 1:size(ICs,1)
    x0 = ICs(i,:);
    [t, x] = ode45(sys, tspan, x0, options);
    plot(x(:,1), x(:,2), 'k', 'LineWidth', 1.5);
    plot(x0(1), x0(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
end

% ---------------------- ͼ������ -------------------------
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

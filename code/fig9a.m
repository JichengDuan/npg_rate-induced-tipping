clc; clear;

% ---------------------- 参数设置 -------------------------
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon = 0.1;
c_start = 0.4;
c_end = 0.1;
Tr = 2;
Tend = 200;
v = (c_end - c_start)/Tr;

% ---------------------- 初值设置 -------------------------
ICs = [
    0, 4;
    0, 2.5;
    1.5, 0;
    1.469 1.469
%    1.5, 0
];

% ---------------------- 网格设置（用于吸引域） -------------------------
x1_vals = linspace(0, 4, 50);
x2_vals = linspace(0, 4, 50);
[X1, X2] = meshgrid(x1_vals, x2_vals);
final_states = nan(size(X1));
attractor_list = [];
tol = 1e-2;

% ---------------------- 吸引域计算 -------------------------
options = odeset('RelTol',1e-6, 'AbsTol',1e-8);
for i = 1:numel(X1)
    x0 = [X1(i); X2(i); c_start];
    [~, sol] = ode45(@(t,x) rate_sys(t,x,s1,s2,s3,epsilon,v,Tr), [0 Tend], x0, options);
    x_final = sol(end,1:2);
    assigned = false;
    for j = 1:size(attractor_list,1)
        if norm(x_final - attractor_list(j,:)) < tol
            final_states(i) = j;
            assigned = true;
            break;
        end
    end
    if ~assigned
        attractor_list = [attractor_list; x_final];
        final_states(i) = size(attractor_list,1);
    end
end

% ---------------------- 图形初始化 -------------------------
figure; hold on;
set(gca,'FontSize',12);
xlabel('$X_1$', 'Interpreter','latex');
ylabel('$X_2$', 'Interpreter','latex');
%title('Combined: Basin, Nullclines, Attractors, and Trajectories', 'Interpreter','latex');
axis([0 4 0 4]); box on; grid on;

% ---------------------- 吸引域背景 -------------------------
imagesc(x1_vals, x2_vals, reshape(final_states, size(X1)), [1 size(attractor_list,1)]);
set(gca,'YDir','normal');
colormap(parula(size(attractor_list,1)));
colorbar;

% ---------------------- nullclines -------------------------
[X1f, X2f] = meshgrid(linspace(0,4,200), linspace(0,4,200));
nullcline1 = - (X1f - s1).*(X1f - s2).*(X1f - s3) + c_end*X2f;
nullcline2 = - epsilon*(X2f - s1).*(X2f - s2).*(X2f - s3);
contour(X1f, X2f, nullcline1, [0 0], 'r', 'LineWidth', 2); % dx1/dt = 0
contour(X1f, X2f, nullcline2, [0 0], 'm', 'LineWidth', 2); % dx2/dt = 0

% ---------------------- 起始与终止 nullcline 轨迹 -------------------------
x1 = linspace(0, 4, 200);
f_start = ((x1 - s1).*(x1 - s2).*(x1 - s3)) / c_start;
f_end   = ((x1 - s1).*(x1 - s2).*(x1 - s3)) / c_end;
plot(x1, f_start, 'b--', 'LineWidth', 1.5);   % 起始 nullcline
%plot(x1, f_end,   'g--', 'LineWidth', 1.5);   % 终态 nullcline

% ---------------------- 平衡点（c=0.1）计算 -------------------------
[x1g, x2g] = meshgrid(linspace(0,4,100), linspace(0,4,100));
attractors = [];
for i = 1:numel(x1g)
    x0 = [x1g(i), x2g(i)];
    [xsol, ~, exitflag] = fsolve(@(x) coupled_system(x,s1,s2,s3,epsilon,c_end), x0, ...
                                 optimoptions('fsolve','Display','off'));
    if exitflag > 0 && all(imag(xsol)==0) && all(xsol>=0 & xsol<=4)
        if isempty(attractors)
            attractors = [attractors; xsol];
        else
            dists = vecnorm(attractors - xsol, 2, 2);
            if all(dists > tol)
                attractors = [attractors; xsol];
            end
        end
    end
end

% ---------------------- 吸引子绘图 -------------------------
for i = 1:size(attractors,1)
    x_star = attractors(i,:);
    J = coupled_jacobian(x_star, s1, s2, s3, epsilon, c_end);
    eigs = eig(J);
    if all(real(eigs) < 0)
        plot(x_star(1), x_star(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 8);
    else
        plot(x_star(1), x_star(2), 'ko', 'MarkerSize', 8);
    end
end

% ---------------------- 动态轨道绘制 -------------------------
arrow_len = 0.15;
for i = 1:size(ICs,1)
    x0 = ICs(i,:);
    [t, x] = ode45(@(t,x) [
        - (x(1)-s1)*(x(1)-s2)*(x(1)-s3) + c_t(t,c_start,c_end,Tr)*x(2);
        - epsilon * (x(2)-s1)*(x(2)-s2)*(x(2)-s3)
    ], [0 Tend], x0, options);

    plot(x(:,1), x(:,2), 'c', 'LineWidth', 1.5);
    plot(x0(1), x0(2), 'ko', 'MarkerFaceColor','c', 'MarkerSize', 6);

    % 添加轨道中段箭头
%     mid_idx = round(length(t)/2);
%     if mid_idx < length(x)
%         dir = x(mid_idx+1,:) - x(mid_idx,:);
%         dir = dir / norm(dir);
%         quiver(x(mid_idx,1), x(mid_idx,2), arrow_len*dir(1), arrow_len*dir(2), ...
%                'MaxHeadSize', 3, 'AutoScale', 'off', ...
%                'Color', 'k', 'LineWidth', 1.2);
%     end
end

% ---------------------- 图例 -------------------------
% legend({'$f_1=0$, $c=0.4$', '$f_1=0$, $c=0.1$', ...
%         '$dX_1/dt=0$', '$dX_2/dt=0$', ...
%         'Stable attractor', 'Unstable point'}, ...
%         'Interpreter','latex', 'Location','northeastoutside');

% ---------------------- 函数定义 -------------------------
function dxdt = rate_sys(t, x, s1, s2, s3, epsilon, v, Tr)
    dxdt = zeros(3,1);
    c = x(3);
    dxdt(1) = - (x(1)-s1)*(x(1)-s2)*(x(1)-s3) + c * x(2);
    dxdt(2) = - epsilon * (x(2)-s1)*(x(2)-s2)*(x(2)-s3);
    dxdt(3) = (t <= Tr) * v;
end

function dxdt = coupled_system(x, s1, s2, s3, epsilon, c)
    dxdt = zeros(2,1);
    dxdt(1) = - (x(1)-s1)*(x(1)-s2)*(x(1)-s3) + c * x(2);
    dxdt(2) = - epsilon * (x(2)-s1)*(x(2)-s2)*(x(2)-s3);
end

function J = coupled_jacobian(x, s1, s2, s3, epsilon, c)
    dfdx1 = - (3*x(1)^2 - 2*(s1+s2+s3)*x(1) + (s1*s2 + s1*s3 + s2*s3));
    dfdx2 = c;
    dgdx1 = 0;
    dgdx2 = - epsilon * (3*x(2)^2 - 2*(s1+s2+s3)*x(2) + (s1*s2 + s1*s3 + s2*s3));
    J = [dfdx1, dfdx2; dgdx1, dgdx2];
end

function c = c_t(t, c0, c1, Tr)
    if t <= Tr
        c = c0 + (c1 - c0) * (t / Tr);
    else
        c = c1;
    end
end

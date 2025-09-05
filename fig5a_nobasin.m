%clc; clear;
hold on
% 参数定义
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon = 1;      % 响应系统时间尺度
c = 0.1283;            % 耦合强度

% 定义状态空间
x1_vals = linspace(0, 4, 100);
x2_vals = linspace(0, 4, 100);
[X1, X2] = meshgrid(x1_vals, x2_vals);

% 零流线计算
nullcline1 = - (X1 - s1).*(X1 - s2).*(X1 - s3) + c*X2;   % dx1/dt = 0
nullcline2 = - epsilon*(X2 - s1).*(X2 - s2).*(X2 - s3);  % dx2/dt = 0

% 零流线绘图
figure(3); hold on;
contour(X1, X2, nullcline1, [0 0], 'r', 'LineWidth', 2); % X1 零流线
contour(X1, X2, nullcline2, [0 0], 'm', 'LineWidth', 2); % X2 零流线

% 自动寻找所有平衡点
[x1g, x2g] = meshgrid(linspace(0,4,100), linspace(0,4,100));
attractors = [];
tol = 1e-2;
for i = 1:numel(x1g)
    x0 = [x1g(i), x2g(i)];
    [xsol, fval, exitflag] = fsolve(@(x) coupled_system(x,s1,s2,s3,epsilon,c), x0, optimoptions('fsolve','Display','off'));
    if exitflag > 0
        if all(imag(xsol)==0) && all(xsol>=0 & xsol<=4)
            is_new = true;
            for j = 1:size(attractors,1)
                if norm(xsol - attractors(j,:)) < tol
                    is_new = false;
                    break;
                end
            end
            if is_new
                attractors = [attractors; xsol];
            end
        end
    end
end

% 稳定性判断并绘制吸引子
for i = 1:size(attractors,1)
    x_star = attractors(i,:);
    J = coupled_jacobian(x_star, s1, s2, s3, epsilon, c);
    eigs = eig(J);
    if all(real(eigs) < 0)
        plot(x_star(1), x_star(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 8); % 稳定吸引子
    else
        plot(x_star(1), x_star(2), 'ko', 'MarkerSize', 8); % 不稳定点
    end
end

xlabel('$X_1$', 'Interpreter','latex'); ylabel('$X_2$', 'Interpreter','latex');
legend({'$dX_1/dt=0$', '$dX_2/dt=0$', 'Stable attractor', 'Unstable fixed point'}, ...
        'Interpreter','latex', 'Location','best');

axis([0 4 0 4]); grid on; box on;


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

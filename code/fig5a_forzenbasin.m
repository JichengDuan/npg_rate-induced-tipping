clc; clear all;

% 参数定义
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon = 1;
c = 0.1;

% --- 吸引域计算 ---
x1_vals = linspace(0, 4, 50);
x2_vals = linspace(0, 4, 50);
[X1, X2] = meshgrid(x1_vals, x2_vals);
final_states = nan(size(X1));
Tend = 200;
options = odeset('RelTol',1e-6, 'AbsTol',1e-8);
attractor_list = [];
tol = 1e-2;

for i = 1:numel(X1)
    x0 = [X1(i), X2(i)];
    [~, sol] = ode45(@(t,x) sys(t,x,s1,s2,s3,epsilon,c), [0 Tend], x0, options);
    x_final = sol(end,:);
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

% --- 吸引域绘图（背景图层） ---
figure(3); clf; hold on;

% 显示吸引域，强制 color range 为 [1, num_attractors]
imagesc(x1_vals, x2_vals, reshape(final_states, size(X1)), [1 size(attractor_list,1)]);
set(gca,'YDir','normal');
colormap(parula(size(attractor_list,1)));
colorbar;

% --- 零流线绘图 ---
x1_vals_dense = linspace(0, 4, 100);
x2_vals_dense = linspace(0, 4, 100);
[X1d, X2d] = meshgrid(x1_vals_dense, x2_vals_dense);
nullcline1 = - (X1d - s1).*(X1d - s2).*(X1d - s3) + c*X2d;
nullcline2 = - epsilon*(X2d - s1).*(X2d - s2).*(X2d - s3);
contour(X1d, X2d, nullcline1, [0 0], 'r', 'LineWidth', 2); % dx1/dt = 0
contour(X1d, X2d, nullcline2, [0 0], 'b', 'LineWidth', 2); % dx2/dt = 0

% --- 平衡点计算并绘图 ---
[x1g, x2g] = meshgrid(linspace(0,4,100), linspace(0,4,100));
attractors = [];
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

for i = 1:size(attractors,1)
    x_star = attractors(i,:);
    J = coupled_jacobian(x_star, s1, s2, s3, epsilon, c);
    eigs = eig(J);
    if all(real(eigs) < 0)
        plot(x_star(1), x_star(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 8); % 稳定
    else
        plot(x_star(1), x_star(2), 'ko', 'MarkerSize', 8); % 不稳定
    end
end

% --- 图形修饰 ---
xlabel('$X_1$', 'Interpreter','latex'); 
ylabel('$X_2$', 'Interpreter','latex');
title(['Combined Phase Portrait & Basins, $c = ', num2str(c), '$'],'Interpreter','latex');
legend({'$dX_1/dt=0$', '$dX_2/dt=0$', 'Stable attractor', 'Unstable fixed point'}, ...
        'Interpreter','latex', 'Location','northeast');
axis([0 4 0 4]); grid on; box on;

% --- 函数定义 ---
function dxdt = sys(t, x, s1, s2, s3, epsilon, c)
    dxdt = zeros(2,1);
    dxdt(1) = - (x(1)-s1)*(x(1)-s2)*(x(1)-s3) + c * x(2);
    dxdt(2) = - epsilon * (x(2)-s1)*(x(2)-s2)*(x(2)-s3);
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

%% Fig. 7(a): Relative size of non-autonomous basins vs Tr (drive-response, Eq. (8))
clear; clc; close all;

% -------------------- Parameters --------------------
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon = 1.0;                 % as in the paper
c_start = 0.4;
c_end   = 0.1;

% 扫描的 Tr（可按需加密/扩展）
Tr_list = linspace(0, 16, 20);     % 覆盖快漂移到慢漂移
Tend_total = 220;                    % 总积分时间（含漂移后收敛）
Tsettle_end = 120;                   % 漂移结束后在 c_end 下的“收敛时间”

% 初值网格（非负约束的状态空间）
Ngrid = 65;                          % 网格密度（可调）
x1_vals = linspace(0, 4, Ngrid);
x2_vals = linspace(0, 4, Ngrid);
[X1g, X2g] = meshgrid(x1_vals, x2_vals);
IC = [X1g(:), X2g(:)];
N_IC = size(IC,1);

% ODE 选项
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

% -------------------- 终态（c_end）下的稳定吸引子 --------------------
% 冻结系统的方程： (X1-s1)(X1-s2)(X1-s3) = c_end * X2
% => X1^3 - S1 X1^2 + S2 X1 - (S3 + c_end*X2) = 0
S1 = s1 + s2 + s3;               % 6
S2 = s1*s2 + s1*s3 + s2*s3;      % 11
S3 = s1*s2*s3;                   % 6
x2_levels = [s1, s2, s3];

% 找到 c_end 下的所有稳定平衡
stable_eqs = [];
for X2eq = x2_levels
    coeffs = [1, -S1, S2, -(S3 + c_end*X2eq)];
    r = roots(coeffs);
    r = r(abs(imag(r)) < 1e-10);  r = real(r);
    for k = 1:numel(r)
        X1eq = r(k);
        % 冻结系统雅可比（上三角），稳定性由对角元决定
        df1dx1 = -(3*X1eq^2 - 2*S1*X1eq + S2);
        df2dx2 = -epsilon*(3*X2eq^2 - 2*S1*X2eq + S2);  % x2=1,3 为负；x2=2 为正(不稳定)
        if (df2dx2 < 0) && (df1dx1 < 0)
            stable_eqs = [stable_eqs; X1eq, X2eq];
        end
    end
end
% 去重并排序（按 X2 升序，再按 X1 升序）以固定 “basin 1..4” 的含义
stable_eqs = unique(round(stable_eqs, 8), 'rows');
[~, idxSort] = sortrows(stable_eqs, [2 1]);
stable_eqs = stable_eqs(idxSort, :);
n_atts = size(stable_eqs,1);
if n_atts ~= 4
    warning('在 c_end=%.3f 下找到了 %d 个稳定平衡（期望为 4）。', c_end, n_atts);
end

% -------------------- 主循环：对每个 Tr 统计相对面积 --------------------
basin_frac = zeros(numel(Tr_list), n_atts);

for k = 1:numel(Tr_list)
    Tr = Tr_list(k);
    v  = (c_end - c_start)/Tr;

    % 非自治系统：c(t) 线性漂移，再保持
    c_t = @(t) (t<=Tr).*(c_start + v*t) + (t>Tr).*c_end;
    rhs_rate = @(t,x) [ ...
        - (x(1)-s1).*(x(1)-s2).*(x(1)-s3) + c_t(t).*x(2); ...
        - epsilon * (x(2)-s1).*(x(2)-s2).*(x(2)-s3) ];

    % 终态（冻结）系统：用于“打磨收敛”
    rhs_end  = @(t,x) [ ...
        - (x(1)-s1).*(x(1)-s2).*(x(1)-s3) + c_end.*x(2); ...
        - epsilon * (x(2)-s1).*(x(2)-s2).*(x(2)-s3) ];

    counts = zeros(1, n_atts);

    for i = 1:N_IC
        x0 = IC(i,:);

        % 先积分非自治段（0->Tend_total）
        [~, X] = ode45(rhs_rate, [0 Tend_total - Tsettle_end], x0, opts);
        x_mid = X(end,:);

        % 再在 c_end 下继续积分 Tsettle_end，确保真正落入某个吸引子盆
        [~, X2] = ode45(rhs_end, [0 Tsettle_end], x_mid, opts);
        xf = X2(end,:);

        % 归类：按“最终收敛到哪个稳定平衡点”
        % （采用最近吸引子规则，但此时已在冻结系统下收敛，鲁棒）
        [~, id_min] = min(sum((stable_eqs - xf).^2, 2));
        counts(id_min) = counts(id_min) + 1;
    end

    basin_frac(k, :) = counts / N_IC;
    fprintf('Tr=%.3f (v=%+.4f):  %s\n', Tr, v, mat2str(basin_frac(k,:),3));
end

% -------------------- 绘图（论文风格的 4 条黑线） --------------------
figure('Color','w'); hold on;
plot(Tr_list, basin_frac(:,1), 'k-',  'LineWidth',1.8);   % basin 1
plot(Tr_list, basin_frac(:,2), 'k--', 'LineWidth',1.8);   % basin 2
plot(Tr_list, basin_frac(:,3), 'k:',  'LineWidth',1.8);   % basin 3
plot(Tr_list, basin_frac(:,4), 'k-.', 'LineWidth',1.8);   % basin 4
xlabel('T_r'); ylabel('relative size of basin');
xlim([min(Tr_list) max(Tr_list)]);
ylim([0 0.52]); grid on; box on;
title('(a) Relative size of non-autonomous basins');
legend({'basin 1','basin 2','basin 3','basin 4'}, 'Location','eastoutside');

%（可选）打印终态吸引子坐标，便于核对
disp('Stable equilibria at c_end (sorted as basin 1..4):');
disp(stable_eqs);

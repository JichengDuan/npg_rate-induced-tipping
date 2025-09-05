%% Fig. 7(b) — tipping probabilities P21 & P43 vs Tr (robust, layer-based)
clear; clc; close all;

% ---------------- Parameters ----------------
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon  = 1.0;
c_start  = 0.4;
c_end    = 0.1;

Tr_list  = linspace(0.2, 16, 22);   % fast -> slow drift

% Numerics
Ngrid    = 81;      % cell-centered grid (odd recommended)
Tfreeze  = 200;     % settle time in frozen systems
TendRate = 220;     % non-autonomous integration time (drift + buffer)
Tsettle  = 120;     % extra settle at c_end
dt0      = 1e-2;    % tiny pre-step
opts     = odeset('RelTol',1e-6,'AbsTol',1e-8);

% ---------------- Cell-centered initial grid (first quadrant) ----------------
x1_edges = linspace(0,4,Ngrid+1);
x2_edges = linspace(0,4,Ngrid+1);
x1c = 0.5*(x1_edges(1:end-1)+x1_edges(2:end));
x2c = 0.5*(x2_edges(1:end-1)+x2_edges(2:end));
[X1c,X2c] = meshgrid(x1c,x2c);
IC = [X1c(:), X2c(:)];
N_IC = size(IC,1);

% ---------------- Helpers ----------------
S1 = s1+s2+s3; S2 = s1*s2 + s1*s3 + s2*s3; S3 = s1*s2*s3;
x2_levels = [s1 s2 s3];

find_stable = @(cfix) local_find_stable(cfix,x2_levels,S1,S2,S3,epsilon);
rhs_frozen  = @(cfix) @(t,x)[ ...
    - (x(1)-s1).*(x(1)-s2).*(x(1)-s3) + cfix.*x(2); ...
    - epsilon*(x(2)-s1).*(x(2)-s2).*(x(2)-s3) ];

% ---------------- Identify start/end equilibria ----------------
EQ_start = find_stable(c_start);   % typically 2 stable points: at X2=1 and X2=3
EQ_end   = find_stable(c_end);     % typically 4 stable points: two at X2=1, two at X2=3

% 起始：找“X2=1 的稳定点”与“X2=3 的稳定点”
idxStart_1 = pick_by_layer(EQ_start, 1.0);    % 唯一
idxStart_3 = pick_by_layer(EQ_start, 3.0);    % 唯一
if isempty(idxStart_1) || isempty(idxStart_3)
    error('在 c_start=%.2f 下未能在 X2=1/3 找到稳定吸引子。', c_start);
end

% 终态：在 X2=1 与 X2=3 层，各找“低 X1（左）”的那个
[idxEnd_1_low,  idxEnd_1_high] = pick_two_by_layer(EQ_end, 1.0);
[idxEnd_3_low,  idxEnd_3_high] = pick_two_by_layer(EQ_end, 3.0);
if isempty(idxEnd_1_low) || isempty(idxEnd_3_low)
    error('在 c_end=%.2f 下未能在 X2=1/3 找到两侧稳定吸引子。', c_end);
end

% ---------------- Step 1: 起始冻结系统下标记“起始层盆” ----------------
start_label = zeros(N_IC,1);   % 1: 来自 X2=1 层盆； 3: 来自 X2=3 层盆； 0: 其它/未识别
for i=1:N_IC
    x = IC(i,:);
    [~,Xp] = ode45(rhs_frozen(c_start), [0 dt0], x, opts); x = Xp(end,:);
    [~,Xs] = ode45(rhs_frozen(c_start), [0 Tfreeze], x, opts);
    xf = Xs(end,:);
    % 归最近的起始稳定点
    [~,idx] = min(sum((EQ_start - xf).^2,2));
    if idx==idxStart_1
        start_label(i) = 1;
    elseif idx==idxStart_3
        start_label(i) = 3;
    else
        start_label(i) = 0;
    end
end
mask1 = (start_label==1);  N1 = sum(mask1);
mask3 = (start_label==3);  N3 = sum(mask3);
fprintf('起始盆计数：X2=1 层 N1=%d,  X2=3 层 N3=%d（总样本 %d）\n', N1, N3, N_IC);

% ---------------- Step 2: 对每个 Tr 计算 P21 与 P43 ----------------
P21 = nan(numel(Tr_list),1);  % from start-layer X2=1  -> end-layer X2=1 low-X1
P43 = nan(numel(Tr_list),1);  % from start-layer X2=3  -> end-layer X2=3 low-X1

for k=1:numel(Tr_list)
    Tr = Tr_list(k);
    v  = (c_end - c_start)/Tr;
    c_t = @(t) (t<=Tr).*(c_start + v*t) + (t>Tr).*c_end;
    rhs_rate = @(t,x)[ ...
        - (x(1)-s1).*(x(1)-s2).*(x(1)-s3) + c_t(t).*x(2); ...
        - epsilon*(x(2)-s1).*(x(2)-s2).*(x(2)-s3) ];

    % 计数
    hit1_to_1low = 0;  % X2=1 -> (end) X2=1 low-X1
    hit3_to_3low = 0;  % X2=3 -> (end) X2=3 low-X1

    for i=1:N_IC
        x = IC(i,:);
        [~,Xp]   = ode45(rhs_frozen(c_start), [0 dt0], x, opts); x = Xp(end,:);
        [~,Xr]   = ode45(rhs_rate, [0 TendRate], x, opts);      x = Xr(end,:);
        [~,Xe]   = ode45(rhs_frozen(c_end), [0 Tsettle], x, opts); xf = Xe(end,:);
        % 终态归类
        [~,idxEnd] = min(sum((EQ_end - xf).^2,2));
        if mask1(i) && (idxEnd==idxEnd_1_low)
            hit1_to_1low = hit1_to_1low + 1;
        elseif mask3(i) && (idxEnd==idxEnd_3_low)
            hit3_to_3low = hit3_to_3low + 1;
        end
    end

    P21(k) = hit1_to_1low / max(N1,1);
    P43(k) = hit3_to_3low / max(N3,1);

    fprintf('Tr=%.3f  P21=%.4f  P43=%.4f\n', Tr, P21(k), P43(k));
end

% ---------------- Plot （与论文风格一致） ----------------
figure('Color','w'); hold on;
plot(Tr_list, P21, 'k-', 'LineWidth', 1.8);   % 实线
plot(Tr_list, P43, 'k:', 'LineWidth', 1.8);   % 点线
xlabel('T_r'); ylabel('tipping probability');
xlim([min(Tr_list) max(Tr_list)]); ylim([0 0.5]);
grid on; box on;
title('(b) tipping probabilities');
legend('P_{21}','P_{43}','Location','northeast');

%% ---------- local functions ----------
function eqs = local_find_stable(cfix,x2_levels,S1,S2,S3,epsilon)
    eqs = [];
    for X2eq = x2_levels
        coeffs = [1, -S1, S2, -(S3 + cfix*X2eq)]; % x1^3 - S1 x1^2 + S2 x1 - (S3 + c X2)=0
        r = roots(coeffs);
        r = r(abs(imag(r))<1e-12); r = real(r);
        for kk = 1:numel(r)
            X1eq = r(kk);
            df1dx1 = -(3*X1eq^2 - 2*S1*X1eq + S2);
            df2dx2 = -epsilon*(3*X2eq^2 - 2*S1*X2eq + S2);
            if df1dx1 < 0 && df2dx2 < 0
                eqs = [eqs; X1eq, X2eq];
            end
        end
    end
    if ~isempty(eqs)
        eqs = unique(round(eqs,8),'rows');
        [~,ord] = sortrows(eqs,[2 1]); % sort by X2 (layer), then X1
        eqs = eqs(ord,:);
    end
end

function idx = pick_by_layer(eqs, x2target)
    tol = 1e-3;
    idxs = find(abs(eqs(:,2)-x2target) < tol);
    if isempty(idxs)
        idx = [];
    else
        % 若有两个，则选更靠右的（通常起始只有一个）
        [~,ord] = sort(eqs(idxs,1),'ascend');
        idx = idxs(ord(end));
    end
end

function [idx_low, idx_high] = pick_two_by_layer(eqs, x2target)
    tol = 1e-3; idx_low=[]; idx_high=[];
    idxs = find(abs(eqs(:,2)-x2target) < tol);
    if numel(idxs)>=2
        [~,ord] = sort(eqs(idxs,1),'ascend');
        idx_low  = idxs(ord(1));   % 左（小X1）
        idx_high = idxs(ord(2));   % 右（大X1）
    end
end

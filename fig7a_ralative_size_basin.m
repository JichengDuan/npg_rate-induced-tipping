%% Fig. 7(a): Relative size of non-autonomous basins vs Tr (drive-response, Eq. (8))
clear; clc; close all;

% -------------------- Parameters --------------------
s1 = 1.0; s2 = 2.0; s3 = 3.0;
epsilon = 1.0;                 % as in the paper
c_start = 0.4;
c_end   = 0.1;

% ɨ��� Tr���ɰ������/��չ��
Tr_list = linspace(0, 16, 20);     % ���ǿ�Ư�Ƶ���Ư��
Tend_total = 220;                    % �ܻ���ʱ�䣨��Ư�ƺ�������
Tsettle_end = 120;                   % Ư�ƽ������� c_end �µġ�����ʱ�䡱

% ��ֵ���񣨷Ǹ�Լ����״̬�ռ䣩
Ngrid = 65;                          % �����ܶȣ��ɵ���
x1_vals = linspace(0, 4, Ngrid);
x2_vals = linspace(0, 4, Ngrid);
[X1g, X2g] = meshgrid(x1_vals, x2_vals);
IC = [X1g(:), X2g(:)];
N_IC = size(IC,1);

% ODE ѡ��
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

% -------------------- ��̬��c_end���µ��ȶ������� --------------------
% ����ϵͳ�ķ��̣� (X1-s1)(X1-s2)(X1-s3) = c_end * X2
% => X1^3 - S1 X1^2 + S2 X1 - (S3 + c_end*X2) = 0
S1 = s1 + s2 + s3;               % 6
S2 = s1*s2 + s1*s3 + s2*s3;      % 11
S3 = s1*s2*s3;                   % 6
x2_levels = [s1, s2, s3];

% �ҵ� c_end �µ������ȶ�ƽ��
stable_eqs = [];
for X2eq = x2_levels
    coeffs = [1, -S1, S2, -(S3 + c_end*X2eq)];
    r = roots(coeffs);
    r = r(abs(imag(r)) < 1e-10);  r = real(r);
    for k = 1:numel(r)
        X1eq = r(k);
        % ����ϵͳ�ſɱȣ������ǣ����ȶ����ɶԽ�Ԫ����
        df1dx1 = -(3*X1eq^2 - 2*S1*X1eq + S2);
        df2dx2 = -epsilon*(3*X2eq^2 - 2*S1*X2eq + S2);  % x2=1,3 Ϊ����x2=2 Ϊ��(���ȶ�)
        if (df2dx2 < 0) && (df1dx1 < 0)
            stable_eqs = [stable_eqs; X1eq, X2eq];
        end
    end
end
% ȥ�ز����򣨰� X2 �����ٰ� X1 �����Թ̶� ��basin 1..4�� �ĺ���
stable_eqs = unique(round(stable_eqs, 8), 'rows');
[~, idxSort] = sortrows(stable_eqs, [2 1]);
stable_eqs = stable_eqs(idxSort, :);
n_atts = size(stable_eqs,1);
if n_atts ~= 4
    warning('�� c_end=%.3f ���ҵ��� %d ���ȶ�ƽ�⣨����Ϊ 4����', c_end, n_atts);
end

% -------------------- ��ѭ������ÿ�� Tr ͳ�������� --------------------
basin_frac = zeros(numel(Tr_list), n_atts);

for k = 1:numel(Tr_list)
    Tr = Tr_list(k);
    v  = (c_end - c_start)/Tr;

    % ������ϵͳ��c(t) ����Ư�ƣ��ٱ���
    c_t = @(t) (t<=Tr).*(c_start + v*t) + (t>Tr).*c_end;
    rhs_rate = @(t,x) [ ...
        - (x(1)-s1).*(x(1)-s2).*(x(1)-s3) + c_t(t).*x(2); ...
        - epsilon * (x(2)-s1).*(x(2)-s2).*(x(2)-s3) ];

    % ��̬�����ᣩϵͳ�����ڡ���ĥ������
    rhs_end  = @(t,x) [ ...
        - (x(1)-s1).*(x(1)-s2).*(x(1)-s3) + c_end.*x(2); ...
        - epsilon * (x(2)-s1).*(x(2)-s2).*(x(2)-s3) ];

    counts = zeros(1, n_atts);

    for i = 1:N_IC
        x0 = IC(i,:);

        % �Ȼ��ַ����ζΣ�0->Tend_total��
        [~, X] = ode45(rhs_rate, [0 Tend_total - Tsettle_end], x0, opts);
        x_mid = X(end,:);

        % ���� c_end �¼������� Tsettle_end��ȷ����������ĳ����������
        [~, X2] = ode45(rhs_end, [0 Tsettle_end], x_mid, opts);
        xf = X2(end,:);

        % ���ࣺ���������������ĸ��ȶ�ƽ��㡱
        % ��������������ӹ��򣬵���ʱ���ڶ���ϵͳ��������³����
        [~, id_min] = min(sum((stable_eqs - xf).^2, 2));
        counts(id_min) = counts(id_min) + 1;
    end

    basin_frac(k, :) = counts / N_IC;
    fprintf('Tr=%.3f (v=%+.4f):  %s\n', Tr, v, mat2str(basin_frac(k,:),3));
end

% -------------------- ��ͼ�����ķ��� 4 �����ߣ� --------------------
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

%����ѡ����ӡ��̬���������꣬���ں˶�
disp('Stable equilibria at c_end (sorted as basin 1..4):');
disp(stable_eqs);

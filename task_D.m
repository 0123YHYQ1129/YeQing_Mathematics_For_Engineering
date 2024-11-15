% 定义参数范围
delta_values = linspace(0, 2, 500);   % delta 从 0 到 1
epsilon_values = linspace(0, 1, 300); % epsilon 从 0 到 0.5
T = 2 * pi;  % 周期 T

% 设置 ODE solver 的选项
options = odeset('AbsTol', 1e-7, 'RelTol', 1e-8);

% 初始化结果矩阵
stability_map = zeros(length(epsilon_values), length(delta_values));

% 定义 A(t) 矩阵
A = @(t, delta, epsilon) [0, 1; -(delta + 2 * epsilon * cos(t)), 0];

% 定义微分方程函数
odefun = @(t, Q_flat, delta, epsilon) ...
    reshape(A(t, delta, epsilon) * reshape(Q_flat, 2, 2), 4, 1);

% 初始条件：单位矩阵的扁平化形式
Q0 = eye(2);
Q0_flat = Q0(:);

% 遍历 epsilon 和 delta 的值
for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    for j = 1:length(delta_values)
        delta = delta_values(j);
        
        % 使用 ode45 进行数值积分，并使用指定的 options
        [~, Q_solution] = ode45(@(t, Q_flat) odefun(t, Q_flat, delta, epsilon), [0 T], Q0_flat, options);
        
        % 提取 Q(T) 并将其转换回 2x2 矩阵形式
        Q_T = reshape(Q_solution(end, :)', 2, 2);
        
        % 计算 Q(T) 的特征值
        eigenvalues = eig(Q_T);
        
        % 判断稳定性
        if max(abs(eigenvalues)) < 1
            stability_map(i, j) = 1; % 稳定
        else
            stability_map(i, j) = 0; % 不稳定
        end
    end
end

% 绘制稳定性图，交换 X 和 Y 轴
figure;
imagesc(delta_values, epsilon_values, stability_map);
set(gca, 'YDir', 'normal');
xlabel('\delta'); % X 轴变为 delta
ylabel('\epsilon'); % Y 轴变为 epsilon
title('Stability Map');
colorbar;
colormap(flipud(gray)); % 使用灰度图，黑色为不稳定，白色为稳定
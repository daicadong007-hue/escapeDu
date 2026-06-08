tic
%% Clear Workspace and Command Window
clc
clear

%% Font size for axes labels and titles
set(groot, 'defaultAxesFontSize', 20);
font_size = 12; % Customize font size here (e.g., 12, 14, 16, etc.)

%% Nhập thông số bánh răng
m = 2e-3;          % Module (m)
N1 = 25;           % Số răng bánh răng nhỏ (pinion)
N2 = 30;           % Số răng bánh răng lớn (gear)
alfa = 20 * pi / 180; % Góc áp lực (rad)
h_a_star = 1;      % Hệ số chiều cao đầu răng
c_star = 0.25;     % Hệ số khe hở
E = 200e9;         % Module Young (Pa)
Nu = 0.3;          % Hệ số Poisson
L = 20e-3;         % Chiều rộng răng (m)
r_int1 = 6.5e-3;   % Bán kính lỗ trục bánh răng nhỏ (m)
r_int2 = 7.6e-3;   % Bán kính lỗ trục bánh răng lớn (m)
r_c_bar= 0.6;      % Bán kính bo tròn đầu dao
%% Nhập thông số vết nứt
q_0_values = [0,0.1e-3, 0.2e-3, 0.3e-3, 0.4e-3, 0.5e-3, 0.6e-3, 0.7e-3, 0.8e-3,0.9e-3, 1e-3, 1.1e-3, 1.2e-3]; % Các độ sâu vết nứt (m)
q_0_healthy = 0;       % Chiều sâu vết nứt cho răng khỏe mạnh (m)
alfa_c = 60 * pi / 180; % Góc nghiêng của vết nứt (rad)
cracked_tooth_index = 1; % Chỉ số của răng bị nứt (bắt đầu từ 1)

%% Thêm hệ số ma sát (theo He 2007)
mu = 0.05; % Hệ số ma sát

%% Thêm hệ số giảm chấn
zeta = 0.07; % Hệ số giảm chấn

%% Tính toán các giá trị ban đầu
G = E / (2 * (1 + Nu)); % Module trượt (Pa)
r1 = m * N1 / 2;       % Bán kính pitch của bánh răng nhỏ (m)
r_a1 = r1 + h_a_star * m; % Bán kính đầu răng
r_f1 = r1 - h_a_star * m; % Bán kính chân răng
r_d1 = r1 - (c_star + h_a_star) * m; % Bán kính đáy răng
r_b1 = r1 * cos(alfa); % Bán kính cơ sở (m)
r2 = m * N2 / 2;       % Bán kính pitch của bánh răng lớn (m)
r_a2 = r2 + h_a_star * m;
r_f2 = r2 - h_a_star * m;
r_d2 = r2 - (c_star + h_a_star) * m;
r_b2 = r2 * cos(alfa);
eps_alfa = (sqrt(r_a2^2 - r_b2^2) + sqrt(r_a1^2 - r_b1^2) - (r1 + r2) * sin(alfa)) / (pi * m * cos(alfa)); % Hệ số trùng lặp
inv_alfa = tan(alfa) - alfa; % Involute của góc áp lực
teta_b1 = pi / (2 * N1) + inv_alfa; % Góc tại cơ sở của bánh răng nhỏ
teta_b2 = pi / (2 * N2) + inv_alfa; % Góc tại cơ sở của bánh răng lớn
h_f1 = r_d1 / r_int1; % Tỷ lệ cho bánh răng nhỏ
h_f2 = r_d2 / r_int2; % Tỷ lệ cho bánh răng lớn

%% Tính các góc ăn khớp
syms x
if r_b1 < r_f1
    alfa_01 = bsm2(sqrt((r_b1 * (x + teta_b1))^2 + r_b1^2) - r_f1, 0, pi/2, 1e-4);
    alfa_12 = bsm2(sqrt((r_b2 * (x + teta_b2))^2 + r_b2^2) - r_a2, 0, pi/2, 1e-4);
else
    alfa_01 = 0;
    alfa_12 = bsm2(sqrt((r_b2 * (x + teta_b2))^2 + r_b2^2) - r_a2 , 0, pi/2, 1e-4);
end
if r_b2 < r_f2
    alfa_02 = bsm2(sqrt((r_b2 * (x + teta_b2))^2 + r_b2^2) - r_f2, 0, pi/2, 1e-4);
    alfa_11 = bsm2(sqrt((r_b1 * (x + teta_b1))^2 + r_b1^2) - r_a1, 0, pi/2, 1e-4);
else
    alfa_02 = 0;
    alfa_11 = bsm2(sqrt((r_b1 * (x + teta_b1))^2 + r_b1^2) - r_a1, 0, pi/2, 1e-4);
end

%% Tính các góc và chiều dài chân răng
if r_b1 < r_d1
    beta_01 = bsm2(r_b1 * ((x + teta_b1) * sin(x) + cos(x)) - r_f1, 0, pi/2, 1e-4);
    L_d1 = r_b1 * ((beta_01 + teta_b1) * cos(beta_01) - sin(beta_01));
else
    beta_01 = 0;
    L_d1 = r_b1 * teta_b1;
end
if r_b2 < r_d2
    beta_02 = bsm2(r_b2 * ((x + teta_b2) * sin(x) + cos(x)) - r_f2, 0, pi/2, 1e-4);
    L_d2 = r_b2 * ((beta_02 + teta_b2) * cos(beta_02) - sin(beta_02));
else
    beta_02 = 0;
    L_d2 = r_b2 * teta_b2;
end
teta_f1 = (1/N1) * (pi/2 + 2 * tan(alfa) * (h_a_star - r_c_bar) + 2 * r_c_bar / cos(alfa)); % Góc chân răng bánh răng nhỏ
teta_f2 = (1/N2) * (pi/2 + 2 * tan(alfa) * (h_a_star - r_c_bar) + 2 * r_c_bar / cos(alfa)); % Góc chân răng bánh răng lớn
S_f1 = 2 * teta_f1 * r_d1;   % Chiều dài chân răng của bánh răng nhỏ
S_f2 = 2 * teta_f2 * r_d2;   % Chiều dài chân răng của bánh răng lớn

%% Tính độ cứng của bánh răng 2 (khỏe mạnh)
if r_d2 <= r_b2
    [K_a2 K_b2 K_s2 K_f2] = toothmesh1(E, L, G, r2, r_b2, r_d2, teta_f2, S_f2, h_f2, teta_b2, alfa_02, alfa_12, beta_02, 0, 0);
else
    [K_a2 K_b2 K_s2 K_f2] = toothmesh2(E, L, G, r2, r_b2, r_d2, teta_f2, S_f2, h_f2, teta_b2, alfa_02, alfa_12, beta_02, 0, 0);
end
K_a2 = fliplr(K_a2);
K_b2 = fliplr(K_b2);
K_s2 = fliplr(K_s2);
K_f2 = fliplr(K_f2);

%% Tính độ cứng của bánh răng 1 khỏe mạnh (không có vết nứt)
if r_d1 <= r_b1
    [K_a1_healthy K_b1_healthy K_s1_healthy K_f1_healthy] = toothmesh1(E, L, G, r1, r_b1, r_d1, teta_f1, S_f1, h_f1, teta_b1, alfa_01, alfa_11, beta_01, q_0_healthy, 0);
else
    [K_a1_healthy K_b1_healthy K_s1_healthy K_f1_healthy] = toothmesh2(E, L, G, r1, r_b1, r_d1, teta_f1, S_f1, h_f1, teta_b1, alfa_01, alfa_11, beta_01, q_0_healthy, 0);
end

%% Tính độ cứng tổng hợp cho trường hợp khỏe mạnh
K_A_healthy = 1 ./ (1 ./ K_a1_healthy + 1 ./ K_a2);
K_B_healthy = 1 ./ (1 ./ K_b1_healthy + 1 ./ K_b2);
K_F_healthy = 1 ./ (1 ./ K_f1_healthy + 1 ./ K_f2);
K_S_healthy = 1 ./ (1 ./ K_s1_healthy + 1 ./ K_s2);
K_h = (pi * E * L) / (4 * (1 - Nu^2));
K_h = ones(1, length(K_a1_healthy)) * K_h;
K_healthy = 1 ./ (1 ./ K_h + 1 ./ K_A_healthy + 1 ./ K_B_healthy + 1 ./ K_S_healthy + 1 ./ K_F_healthy);

PTH = length(K_healthy);
K2_healthy = [zeros(1, PTH) K_healthy zeros(1, PTH)];
Pb = floor(PTH / eps_alfa);

%% Tính toán các điểm trong chu kỳ ăn khớp (theo He 2007)
lambda = pi * m * cos(alfa); % Base pitch
omega_p = 2400 * 2 * pi / 60; % Tốc độ góc pinion (rad/s)
t_c = lambda / (omega_p * r_b1); % Thời gian chu kỳ ăn khớp
L_AP = r_b1 * (tan(acos(r_b2 / r_a2)) - tan(alfa)); % Khoảng cách đến pitch point
L_AB = eps_alfa * lambda; % Khoảng cách ăn khớp đôi
t_P = (L_AP / lambda) * t_c; % Thời điểm pitch point
t_B = (L_AB / lambda) * t_c; % Thời điểm kết thúc ăn khớp đôi
pos_double = t_B / t_c; % Tỷ lệ ăn khớp đôi
pos_single = 1; % Tỷ lệ ăn khớp đơn

%% Thông số cho mô hình động lực học
J_p = 0.96e-4;        % Moment quán tính của bánh răng nhỏ (kg·m²)
J_g = 2e-4;           % Moment quán tính của bánh răng lớn (kg·m²)
m_p = 0.3083;         % Khối lượng của bánh răng nhỏ (kg)
m_g = 0.4439;         % Khối lượng của bánh răng lớn (kg)
K_px = 6.56e9;        % Độ cứng ổ trục theo x (N/m)
K_gx = 6.56e9;
K_py = 6.56e9;        % Độ cứng ổ trục theo y (N/m)
K_gy = 6.56e9;
C_px = 1.8e3;         % Damping ổ trục theo x (Ns/m)
C_gx = 1.8e3;
C_py = 1.8e3;         % Damping ổ trục theo y (Ns/m)
C_gy = 1.8e3;
T_g = 60;             % Mô-men xoắn trên bánh răng lớn (Nm)
T_p = T_g * (N1 / N2); % Mô-men xoắn trên bánh răng nhỏ (Nm)

%% Tính moment quán tính hiệu quả
J_e = J_p * J_g / (J_p * r_b2^2 + J_g * r_b1^2);

%% Cấu trúc tham số cho ODE
params.J_p = J_p;
params.J_g = J_g;
params.m_p = m_p;
params.m_g = m_g;
params.r_b1 = r_b1;
params.r_b2 = r_b2;
params.K_px = K_px;
params.K_gx = K_gx;
params.K_py = K_py;
params.K_gy = K_gy;
params.C_px = C_px;
params.C_gx = C_gx;
params.C_py = C_py;
params.C_gy = C_gy;
params.T_p = T_p;
params.T_g = T_g;
params.mu = mu;
params.zeta = zeta;
params.J_e = J_e;
params.t_c = t_c;
params.t_P = t_P;
params.lambda = lambda;
params.theta_period = 360 / N1;
params.N1 = N1;
params.K_healthy = K_healthy;
params.PTH = PTH;
params.overlap_start_index = Pb;
params.cracked_tooth_index = cracked_tooth_index;

%% Điều kiện ban đầu và khoảng thời gian mô phỏng
z0 = [0; 0; omega_p; omega_p * (N1 / N2); 0; 0; 0; 0; 0; 0; 0; 0];
tspan = [0 0.5]; % Mô phỏng trong 0.5 giây

%% Mô phỏng và tính toán cho từng độ sâu vết nứt
num_cases = length(q_0_values);
y_p_data = cell(num_cases, 1);
t_filtered_data = cell(num_cases, 1);
kurtosis_values = zeros(num_cases, 1);
rms_values = zeros(num_cases, 1);

for i = 1:num_cases
    % Tính độ cứng của bánh răng 1 bị nứt cho độ sâu hiện tại
    q_0_cracked = q_0_values(i);
    if r_d1 <= r_b1
        [K_a1_cracked K_b1_cracked K_s1_cracked K_f1_cracked] = toothmesh1(E, L, G, r1, r_b1, r_d1, teta_f1, S_f1, h_f1, teta_b1, alfa_01, alfa_11, beta_01, q_0_cracked, alfa_c);
    else
        [K_a1_cracked K_b1_cracked K_s1_cracked K_f1_cracked] = toothmesh2(E, L, G, r1, r_b1, r_d1, teta_f1, S_f1, h_f1, teta_b1, alfa_01, alfa_11, beta_01, q_0_cracked, alfa_c);
    end
    
    % Tính độ cứng tổng hợp cho trường hợp răng bị nứt
    K_A_cracked = 1 ./ (1 ./ K_a1_cracked + 1 ./ K_a2);
    K_B_cracked = 1 ./ (1 ./ K_b1_cracked + 1 ./ K_b2);
    K_F_cracked = 1 ./ (1 ./ K_f1_cracked + 1 ./ K_f2);
    K_S_cracked = 1 ./ (1 ./ K_s1_cracked + 1 ./ K_s2);
    K_cracked = 1 ./ (1 ./ K_h + 1 ./ K_A_cracked + 1 ./ K_B_cracked + 1 ./ K_S_cracked + 1 ./ K_F_cracked);
    
    % Cập nhật params với K_cracked
    params.K_cracked = K_cracked;
    
    % Giải ODE
    [t, z] = ode45(@(t, z) gear_dyn(t, z, params), tspan, z0);
    y_p = z(:,7); % Dịch chuyển y của bánh răng nhỏ
    
    % Lọc dữ liệu từ 0.2 đến 0.40 giây
    idx = (t >= 0.2) & (t <= 0.40);
    t_filtered = t(idx);
    y_p_filtered = y_p(idx);
    
    % Loại bỏ thành phần DC
    y_p_filtered = y_p_filtered - mean(y_p_filtered);
    
    % Lưu dữ liệu
    y_p_data{i} = y_p_filtered;
    t_filtered_data{i} = t_filtered;
    
    % Tính kurtosis và RMS
    kurtosis_values(i) = calculate_kurtosis(y_p_filtered);
    rms_values(i) = calculate_rms(y_p_filtered);
end

%% Tính phần trăm thay đổi so với trường hợp khỏe mạnh
kurtosis_healthy = kurtosis_values(1); % Giá trị kurtosis của trường hợp khỏe mạnh (q_0 = 0)
rms_healthy = rms_values(1); % Giá trị RMS của trường hợp khỏe mạnh (q_0 = 0)
kurtosis_percent_change = ((kurtosis_values - kurtosis_healthy) / kurtosis_healthy) * 100;
rms_percent_change = ((rms_values - rms_healthy) / rms_healthy) * 100;

%% Vẽ đồ thị chuyển vị y và phổ tần số riêng rẽ cho từng trường hợp
for i = 1:num_cases
    % Lấy dữ liệu
    t_filtered = t_filtered_data{i};
    y_p_filtered = y_p_data{i};
    
    % Tính toán phổ tần số
    dt = mean(diff(t_filtered));
    Fs = 1 / dt;
    N = length(y_p_filtered);
    f = (0:N-1)*(Fs/N);
    Y = fft(y_p_filtered);
    amplitude = abs(Y) / N;
    amplitude = amplitude(1:floor(N/2)+1);
    f_plot = f(1:floor(N/2)+1);
    
    % Tạo figure cho đồ thị chuyển vị y
    figure('Name', sprintf('Y Displacement for Crack Depth %.1f mm', q_0_values(i)*1e3));
    plot(t_filtered, y_p_filtered * 1e6, 'b-', 'LineWidth', 1.5);
    title(sprintf('Y Displacement (q_0 = %.1f mm)', q_0_values(i)*1e3), 'FontSize', 20);
    xlabel('Time (s)', 'FontSize', 20);
    ylabel('Displacement (μm)', 'FontSize', 20);
    xlim([0.2 0.30]);
    grid on;
    
    % Tạo figure cho đồ thị phổ tần số
    figure('Name', sprintf('Frequency Spectrum for Crack Depth %.1f mm', q_0_values(i)*1e3));
    plot(f_plot, amplitude, 'r-', 'LineWidth', 1.5);
    title(sprintf('Frequency Spectrum (q_0 = %.1f mm)', q_0_values(i)*1e3), 'FontSize', 20);
    xlabel('Frequency (Hz)', 'FontSize', 30);
    ylabel('Amplitude', 'FontSize', 30);
    xlim([0 10000]);
    grid on;
end

%% Vẽ đồ thị so sánh phần trăm thay đổi của Kurtosis và RMS
figure('Name', 'Percentage Change in Kurtosis and RMS vs. Crack Depth');
plot(q_0_values * 1e3, kurtosis_percent_change, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'Kurtosis Change (%)');
hold on;
plot(q_0_values * 1e3, rms_percent_change, 'rs-', 'LineWidth', 1.5, 'DisplayName', 'RMS Change (%)');
title('Percentage Change in Kurtosis and RMS vs. Crack Depth', 'FontSize', font_size);
xlabel('Crack Depth (mm)', 'FontSize', 30);
ylabel('Percentage Change (%)', 'FontSize', 30);
grid on;
legend('show', 'FontSize', 30);
hold off;

toc

%% Hàm tính kurtosis
function kurt = calculate_kurtosis(data)
    n = length(data);
    mean_data = mean(data);
    std_data = std(data);
    kurt = (sum((data - mean_data).^4) / n) / std_data^4;
end

%% Hàm tính RMS
function rms_value = calculate_rms(data)
    x_bar = mean(data);
    rms_value = sqrt(mean((data-x_bar).^2));
end

%% Hàm ODE cho động lực học với một răng bị nứt
function dz = gear_dyn(t, z, params)
    J_p = params.J_p;
    J_g = params.J_g;
    m_p = params.m_p;
    m_g = params.m_g;
    r_b1 = params.r_b1;
    r_b2 = params.r_b2;
    K_px = params.K_px;
    K_gx = params.K_gx;
    K_py = params.K_py;
    K_gy = params.K_gy;
    C_px = params.C_px;
    C_gx = params.C_gx;
    C_py = params.C_py;
    C_gy = params.C_gy;
    T_p = params.T_p;
    T_g = params.T_g;
    mu = params.mu;
    zeta = params.zeta;
    J_e = params.J_e;
    t_c = params.t_c;
    t_P = params.t_P;
    lambda = params.lambda;
    theta_period = params.theta_period;
    N1 = params.N1;
    K_healthy = params.K_healthy;
    K_cracked = params.K_cracked;
    PTH = params.PTH;
    overlap_start_index = params.overlap_start_index;
    cracked_tooth_index = params.cracked_tooth_index;

    theta_p = z(1);
    theta_g = z(2);
    dot_theta_p = z(3);
    dot_theta_g = z(4);
    x_p = z(5);
    x_g = z(6);
    y_p = z(7);
    y_g = z(8);
    dot_x_p = z(9);
    dot_x_g = z(10);
    dot_y_p = z(11);
    dot_y_g = z(12);

    delta = r_b1 * theta_p - r_b2 * theta_g + (y_p - y_g);
    dot_delta = r_b1 * dot_theta_p - r_b2 * dot_theta_g + (dot_y_p - dot_y_g);

    % Xác định răng đang ăn khớp dựa trên góc quay
    theta_p_deg = mod(theta_p * 180 / pi, theta_period * N1);
    tooth_index = floor(theta_p_deg / theta_period) + 1;

    % Điều chỉnh độ cứng dựa trên răng bị nứt
    idx = mod(theta_p_deg / theta_period, 1) * (PTH - 1) + 1;
    if tooth_index == cracked_tooth_index
        k_primary = interp1(1:PTH, K_cracked, idx, 'linear', 'extrap');
    else
        k_primary = interp1(1:PTH, K_healthy, idx, 'linear', 'extrap');
    end
    % Tính damping riêng cho cặp răng chính
    c_primary = 2 * zeta * sqrt(k_primary * J_e);

    % Lực pháp tuyến cho cặp răng chính
    N_primary = max(k_primary * delta + c_primary * dot_delta, 0);

    % Chỉ thêm cặp răng liền kề trong vùng ăn khớp đôi theo hệ số trùng khớp
    N_adjacent = 0;
    if idx > overlap_start_index
        adjacent_idx = idx - overlap_start_index;
        k_adjacent = interp1(1:PTH, K_healthy, adjacent_idx, 'linear', 'extrap'); % Răng liền kề khỏe mạnh
        c_adjacent = 2 * zeta * sqrt(k_adjacent * J_e);
        N_adjacent = max(k_adjacent * delta + c_adjacent * dot_delta, 0);
    end

    % Lực ma sát (dấu dựa trên pitch point)
    t_mod = mod(t, t_c);
    sign_friction = sign(t_mod - t_P); % Đảo chiều tại pitch point
    F_f_primary = mu * N_primary * sign_friction;
    F_f_adjacent = mu * N_adjacent * sign_friction;
    F_f_total = F_f_primary + F_f_adjacent;
    
    % Mô-men do lực pháp tuyến và ma sát
    M_pN_total = (N_primary + N_adjacent) * r_b1;
    M_gN_total = (N_primary + N_adjacent) * r_b2;
    M_pf_total = F_f_total * r_b1;
    M_gf_total = -F_f_total * r_b2;

    dz = zeros(12, 1);
    dz(1) = dot_theta_p;
    dz(2) = dot_theta_g;
    dz(3) = (T_p - M_pN_total + M_pf_total) / J_p;
    dz(4) = (-T_g + M_gN_total + M_gf_total) / J_g;
    dz(5) = dot_x_p;
    dz(6) = dot_x_g;
    dz(7) = dot_y_p;
    dz(8) = dot_y_g;
    dz(9) = (-K_px * x_p - C_px * dot_x_p + F_f_total) / m_p;
    dz(10) = (-K_gx * x_g - C_gx * dot_x_g + F_f_total) / m_g;
    dz(11) = (-(N_primary + N_adjacent) - K_py * y_p - C_py * dot_y_p) / m_p;
    dz(12) = ((N_primary + N_adjacent) - K_gy * y_g - C_gy * dot_y_g) / m_g;
end

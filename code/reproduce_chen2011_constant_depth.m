tic
clc
clear
close all

base_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(base_dir, '..', 'outputs', 'chen2011_constant_depth');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

params = chen2011_parameters();
geom = chen2011_geometry(params);

fprintf('Chen 2011 constant through-width crack reproduction\n');
fprintf('Contact ratio: %.4f\n', geom.eps_alpha);
fprintf('Mesh frequency: %.2f Hz\n', params.omega_p / (2*pi) * params.N_p);

healthy = compute_pair_stiffness(params, geom, 0, params.alpha_c);
crack_depths = (0:0.3:1.2) * 1e-3;
num_cases = numel(crack_depths);

single_profiles = zeros(num_cases, geom.n_mesh);
total_profiles = zeros(num_cases, geom.n_mesh);
time_mesh_profiles = zeros(num_cases, geom.n_rotation);
response_cases = struct([]);

for i = 1:num_cases
    q0 = crack_depths(i);
    cracked = compute_pair_stiffness(params, geom, q0, params.alpha_c);
    single_profiles(i, :) = cracked.K_single;
    total_profiles(i, :) = assemble_total_mesh_profile(healthy.K_single, cracked.K_single, params.cracked_tooth_index, params.cracked_tooth_index, params.N_p, geom);
    time_mesh_profiles(i, :) = rotational_period_mesh_profile(params, geom, healthy.K_single, cracked.K_single);
end

validation = run_table3_validation(params, geom);
disp(validation)

tspan = [0, 0.5];
z0 = [0; 0; params.omega_p; params.omega_g; 0; 0; 0; 0; 0; 0; 0; 0];
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', geom.t_c / 25);

rms_values = zeros(num_cases, 1);
kurtosis_values = zeros(num_cases, 1);
residual_rms_values = zeros(num_cases, 1);
residual_kurtosis_values = nan(num_cases, 1);
mesh_sideband_amp = zeros(num_cases, 1);
mesh_harmonic_amp = zeros(num_cases, 1);

for i = 1:num_cases
    q0 = crack_depths(i);
    cracked = compute_pair_stiffness(params, geom, q0, params.alpha_c);

    dyn = params;
    dyn.geom = geom;
    dyn.K_healthy_single = healthy.K_single;
    dyn.K_cracked_single = cracked.K_single;

    fprintf('Solving ODE for q0 = %.1f mm...\n', q0 * 1e3);
    [t, z] = ode45(@(t, z) gear_rhs(t, z, dyn), tspan, z0, ode_opts);

    keep = t >= 0.22 & t <= 0.40;
    t_ss = t(keep);
    y_p = z(keep, 7);
    y_p = y_p - mean(y_p);

    [freq, amp] = spectrum_from_signal(t_ss, y_p);
    rms_values(i) = sqrt(mean(y_p.^2));
    kurtosis_values(i) = mean(y_p.^4) / mean(y_p.^2)^2;

    mesh_frequency = params.omega_p / (2*pi) * params.N_p;
    [mesh_harmonic_amp(i), mesh_sideband_amp(i)] = spectrum_markers(freq, amp, mesh_frequency, params.omega_p / (2*pi));

    response_cases(i).q0_mm = q0 * 1e3;
    response_cases(i).t = t_ss;
    response_cases(i).y_p = y_p;
    response_cases(i).freq = freq;
    response_cases(i).amp = amp;
end

healthy_time = response_cases(1).t;
healthy_displacement = response_cases(1).y_p;
for i = 2:num_cases
    healthy_on_case_grid = interp1(healthy_time, healthy_displacement, response_cases(i).t, 'linear', 'extrap');
    residual = response_cases(i).y_p - healthy_on_case_grid;
    residual_rms_values(i) = sqrt(mean(residual.^2));
    residual_kurtosis_values(i) = mean(residual.^4) / mean(residual.^2)^2;
end

rms_change = percent_change_from_healthy(rms_values);
kurtosis_change = percent_change_from_healthy(kurtosis_values);

save(fullfile(out_dir, 'chen2011_constant_depth_results.mat'), ...
    'params', 'geom', 'crack_depths', 'single_profiles', 'total_profiles', 'time_mesh_profiles', ...
    'response_cases', 'rms_values', 'kurtosis_values', 'rms_change', ...
    'kurtosis_change', 'residual_rms_values', 'residual_kurtosis_values', ...
    'mesh_harmonic_amp', 'mesh_sideband_amp', 'validation');

writetable(validation, fullfile(out_dir, 'table3_validation.csv'));
write_summary_csv(out_dir, crack_depths, rms_values, kurtosis_values, rms_change, kurtosis_change, residual_rms_values, residual_kurtosis_values, mesh_harmonic_amp, mesh_sideband_amp);
write_time_mesh_csv(out_dir, params, geom, crack_depths, time_mesh_profiles);
make_figures(out_dir, params, geom, crack_depths, single_profiles, total_profiles, time_mesh_profiles, response_cases, rms_change, kurtosis_change, residual_rms_values, residual_kurtosis_values, mesh_harmonic_amp, mesh_sideband_amp);

fprintf('Saved results to: %s\n', out_dir);
toc

function params = chen2011_parameters()
params.m = 2e-3;
params.N_p = 30;
params.N_g = 25;
params.alpha = deg2rad(20);
params.h_a_star = 1;
params.c_star = 0.25;
params.E = 2e11;
params.nu = 0.3;
params.G = params.E / (2 * (1 + params.nu));
params.W = 20e-3;
params.r_int_p = 7.6e-3;
params.r_int_g = 6.5e-3;
params.r_c_bar = 0.6;
params.alpha_c = deg2rad(60);
params.cracked_tooth_index = 1;

params.omega_p = 2000 * 2*pi / 60;
params.omega_g = params.omega_p * params.N_p / params.N_g;
params.T_g = 60;
params.T_p = params.T_g * params.N_p / params.N_g;
params.J_p = 2e-3;
params.J_g = 0.96e-4;
params.m_p = 0.4439;
params.m_g = 0.3083;
params.K_px = 6.56e8;
params.K_gx = 6.56e8;
params.K_py = 6.56e8;
params.K_gy = 6.56e8;
params.C_px = 1.8e3;
params.C_gx = 1.8e3;
params.C_py = 1.8e3;
params.C_gy = 1.8e3;
params.C_m = 67;
params.mu = 0.06;
end

function geom = chen2011_geometry(params)
geom.n_mesh = 1001;

geom.r_p = params.m * params.N_p / 2;
geom.r_ap = geom.r_p + params.h_a_star * params.m;
geom.r_fp = geom.r_p - params.h_a_star * params.m;
geom.r_dp = geom.r_p - (params.c_star + params.h_a_star) * params.m;
geom.r_bp = geom.r_p * cos(params.alpha);

geom.r_g = params.m * params.N_g / 2;
geom.r_ag = geom.r_g + params.h_a_star * params.m;
geom.r_fg = geom.r_g - params.h_a_star * params.m;
geom.r_dg = geom.r_g - (params.c_star + params.h_a_star) * params.m;
geom.r_bg = geom.r_g * cos(params.alpha);

geom.eps_alpha = (sqrt(geom.r_ag^2 - geom.r_bg^2) + sqrt(geom.r_ap^2 - geom.r_bp^2) - ...
    (geom.r_p + geom.r_g) * sin(params.alpha)) / (pi * params.m * cos(params.alpha));

inv_alpha = tan(params.alpha) - params.alpha;
geom.theta_bp = pi / (2 * params.N_p) + inv_alpha;
geom.theta_bg = pi / (2 * params.N_g) + inv_alpha;

geom.h_fp = geom.r_dp / params.r_int_p;
geom.h_fg = geom.r_dg / params.r_int_g;

if geom.r_bp < geom.r_fp
    geom.alpha_0p = bisect_scalar(@(x) sqrt((geom.r_bp * (x + geom.theta_bp))^2 + geom.r_bp^2) - geom.r_fp, 0, pi/2, 1e-10);
else
    geom.alpha_0p = 0;
end
geom.alpha_1p = bisect_scalar(@(x) sqrt((geom.r_bp * (x + geom.theta_bp))^2 + geom.r_bp^2) - geom.r_ap, 0, pi/2, 1e-10);

if geom.r_bg < geom.r_fg
    geom.alpha_0g = bisect_scalar(@(x) sqrt((geom.r_bg * (x + geom.theta_bg))^2 + geom.r_bg^2) - geom.r_fg, 0, pi/2, 1e-10);
else
    geom.alpha_0g = 0;
end
geom.alpha_1g = bisect_scalar(@(x) sqrt((geom.r_bg * (x + geom.theta_bg))^2 + geom.r_bg^2) - geom.r_ag, 0, pi/2, 1e-10);

if geom.r_bp < geom.r_dp
    geom.beta_0p = bisect_scalar(@(x) geom.r_bp * ((x + geom.theta_bp) * sin(x) + cos(x)) - geom.r_fp, 0, pi/2, 1e-10);
else
    geom.beta_0p = 0;
end

if geom.r_bg < geom.r_dg
    geom.beta_0g = bisect_scalar(@(x) geom.r_bg * ((x + geom.theta_bg) * sin(x) + cos(x)) - geom.r_fg, 0, pi/2, 1e-10);
else
    geom.beta_0g = 0;
end

geom.theta_fp = (1 / params.N_p) * (pi/2 + 2 * tan(params.alpha) * (params.h_a_star - params.r_c_bar) + 2 * params.r_c_bar / cos(params.alpha));
geom.theta_fg = (1 / params.N_g) * (pi/2 + 2 * tan(params.alpha) * (params.h_a_star - params.r_c_bar) + 2 * params.r_c_bar / cos(params.alpha));
geom.S_fp = 2 * geom.theta_fp * geom.r_dp;
geom.S_fg = 2 * geom.theta_fg * geom.r_dg;

geom.K_h = pi * params.E * params.W / (4 * (1 - params.nu^2));
geom.base_pitch = pi * params.m * cos(params.alpha);
geom.t_c = geom.base_pitch / (params.omega_p * geom.r_bp);
geom.t_rotation = 2*pi / params.omega_p;
geom.tooth_period = 2*pi / params.N_p;
geom.overlap_start_index = floor(geom.n_mesh / geom.eps_alpha);
geom.n_rotation = params.N_p * (geom.n_mesh - 1) + 1;

geom.L_AP = geom.r_bp * (tan(acos(geom.r_bg / geom.r_ag)) - tan(params.alpha));
geom.t_P = (geom.L_AP / geom.base_pitch) * geom.t_c;
end

function result = compute_pair_stiffness(params, geom, q0, alpha_c)
if geom.r_dp <= geom.r_bp
    [K_ap, K_bp, K_sp, K_fp] = toothmesh1(params.E, params.W, params.G, geom.r_p, geom.r_bp, geom.r_dp, ...
        geom.theta_fp, geom.S_fp, geom.h_fp, geom.theta_bp, geom.alpha_0p, geom.alpha_1p, geom.beta_0p, q0, alpha_c);
else
    [K_ap, K_bp, K_sp, K_fp] = toothmesh2(params.E, params.W, params.G, geom.r_p, geom.r_bp, geom.r_dp, ...
        geom.theta_fp, geom.S_fp, geom.h_fp, geom.theta_bp, geom.alpha_0p, geom.alpha_1p, geom.beta_0p, q0, alpha_c);
end

if geom.r_dg <= geom.r_bg
    [K_ag, K_bg, K_sg, K_fg] = toothmesh1(params.E, params.W, params.G, geom.r_g, geom.r_bg, geom.r_dg, ...
        geom.theta_fg, geom.S_fg, geom.h_fg, geom.theta_bg, geom.alpha_0g, geom.alpha_1g, geom.beta_0g, 0, 0);
else
    [K_ag, K_bg, K_sg, K_fg] = toothmesh2(params.E, params.W, params.G, geom.r_g, geom.r_bg, geom.r_dg, ...
        geom.theta_fg, geom.S_fg, geom.h_fg, geom.theta_bg, geom.alpha_0g, geom.alpha_1g, geom.beta_0g, 0, 0);
end

K_ag = fliplr(K_ag);
K_bg = fliplr(K_bg);
K_sg = fliplr(K_sg);
K_fg = fliplr(K_fg);

K_A = 1 ./ (1 ./ K_ap + 1 ./ K_ag);
K_B = 1 ./ (1 ./ K_bp + 1 ./ K_bg);
K_S = 1 ./ (1 ./ K_sp + 1 ./ K_sg);
K_F = 1 ./ (1 ./ K_fp + 1 ./ K_fg);
K_H = ones(1, numel(K_A)) * geom.K_h;

result.K_single = 1 ./ (1 ./ K_H + 1 ./ K_A + 1 ./ K_B + 1 ./ K_S + 1 ./ K_F);
result.K_A = K_A;
result.K_B = K_B;
result.K_S = K_S;
result.K_F = K_F;
end

function K_total = assemble_total_mesh_profile(K_healthy, K_cracked, cracked_primary_tooth, active_tooth, N_p, geom)
K_total = zeros(1, geom.n_mesh);
for i = 1:geom.n_mesh
    if active_tooth == cracked_primary_tooth
        K_total(i) = K_cracked(i);
    else
        K_total(i) = K_healthy(i);
    end

    if i > geom.overlap_start_index
        adjacent_tooth = mod(active_tooth, N_p) + 1;
        adjacent_idx = i - geom.overlap_start_index;
        if adjacent_tooth == cracked_primary_tooth
            K_total(i) = K_total(i) + K_cracked(adjacent_idx);
        else
            K_total(i) = K_total(i) + K_healthy(adjacent_idx);
        end
    end
end
end

function K_time = rotational_period_mesh_profile(params, geom, K_healthy, K_cracked)
dyn = params;
dyn.geom = geom;
dyn.K_healthy_single = K_healthy;
dyn.K_cracked_single = K_cracked;

theta_samples = linspace(0, 2*pi, geom.n_rotation);
K_time = zeros(1, geom.n_rotation);
for i = 1:geom.n_rotation
    K_time(i) = mesh_stiffness_at_theta(theta_samples(i), dyn);
end
end

function dz = gear_rhs(t, z, params)
geom = params.geom;

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

delta = geom.r_bp * theta_p - geom.r_bg * theta_g + y_p - y_g;
dot_delta = geom.r_bp * dot_theta_p - geom.r_bg * dot_theta_g + dot_y_p - dot_y_g;

k_mesh = mesh_stiffness_at_theta(theta_p, params);
N = max(k_mesh * delta + params.C_m * dot_delta, 0);

t_mod = mod(t, geom.t_c);
friction_sign = sign(t_mod - geom.t_P);
F_f = params.mu * N * friction_sign;

dz = zeros(12, 1);
dz(1) = dot_theta_p;
dz(2) = dot_theta_g;
dz(3) = (params.T_p - N * geom.r_bp + F_f * geom.r_bp) / params.J_p;
dz(4) = (-params.T_g + N * geom.r_bg - F_f * geom.r_bg) / params.J_g;
dz(5) = dot_x_p;
dz(6) = dot_x_g;
dz(7) = dot_y_p;
dz(8) = dot_y_g;
dz(9) = (F_f - params.K_px * x_p - params.C_px * dot_x_p) / params.m_p;
dz(10) = (-F_f - params.K_gx * x_g - params.C_gx * dot_x_g) / params.m_g;
dz(11) = (-N - params.K_py * y_p - params.C_py * dot_y_p) / params.m_p;
dz(12) = (N - params.K_gy * y_g - params.C_gy * dot_y_g) / params.m_g;
end

function k_mesh = mesh_stiffness_at_theta(theta_p, params)
geom = params.geom;
tooth_progress = mod(theta_p / geom.tooth_period, params.N_p);
active_tooth = floor(tooth_progress) + 1;
phase = tooth_progress - floor(tooth_progress);
idx = phase * (geom.n_mesh - 1) + 1;

if active_tooth == params.cracked_tooth_index
    k_mesh = interp1(1:geom.n_mesh, params.K_cracked_single, idx, 'linear');
else
    k_mesh = interp1(1:geom.n_mesh, params.K_healthy_single, idx, 'linear');
end

if idx > geom.overlap_start_index
    adjacent_tooth = mod(active_tooth, params.N_p) + 1;
    adjacent_idx = idx - geom.overlap_start_index;
    if adjacent_tooth == params.cracked_tooth_index
        k_mesh = k_mesh + interp1(1:geom.n_mesh, params.K_cracked_single, adjacent_idx, 'linear');
    else
        k_mesh = k_mesh + interp1(1:geom.n_mesh, params.K_healthy_single, adjacent_idx, 'linear');
    end
end
end

function validation = run_table3_validation(params, geom)
cases = [0, 60, 1.52e8; 0.3e-3, 33, 1.47e8; 0.66e-3, 70, 1.38e8];
labels = ["Healthy"; "Crack No. 1"; "Crack No. 2"];
computed_min = zeros(size(cases, 1), 1);
computed_mean = zeros(size(cases, 1), 1);
computed_max = zeros(size(cases, 1), 1);
target = cases(:, 3);
for i = 1:size(cases, 1)
    stiffness = compute_pair_stiffness(params, geom, cases(i, 1), deg2rad(cases(i, 2)));
    computed_min(i) = min(stiffness.K_single);
    computed_mean(i) = mean(stiffness.K_single);
    computed_max(i) = max(stiffness.K_single);
end
error_pct = 100 * (computed_min - target) ./ target;
validation = table(labels, computed_min, computed_mean, computed_max, target, error_pct, ...
    'VariableNames', {'Case', 'ComputedMin_N_per_m', 'ComputedMean_N_per_m', 'ComputedMax_N_per_m', 'Paper_N_per_m', 'MinError_pct'});
end

function [freq, amp] = spectrum_from_signal(t, x)
dt = mean(diff(t));
Fs = 1 / dt;
N = numel(x);
Y = fft(x);
amp = abs(Y) / N;
amp = amp(1:floor(N/2)+1);
freq = (0:floor(N/2)) * Fs / N;
end

function [harmonic_amp, sideband_amp] = spectrum_markers(freq, amp, mesh_frequency, shaft_frequency)
target_harmonic = 4 * mesh_frequency;
target_sideband = target_harmonic - shaft_frequency;
harmonic_amp = nearest_amp(freq, amp, target_harmonic);
sideband_amp = nearest_amp(freq, amp, target_sideband);
end

function value = nearest_amp(freq, amp, target)
[~, idx] = min(abs(freq - target));
value = amp(idx);
end

function out = percent_change_from_healthy(values)
baseline = values(1);
out = 100 * (values - baseline) / baseline;
end

function write_summary_csv(out_dir, crack_depths, rms_values, kurtosis_values, rms_change, kurtosis_change, residual_rms_values, residual_kurtosis_values, mesh_harmonic_amp, mesh_sideband_amp)
summary = table(crack_depths(:) * 1e3, rms_values, kurtosis_values, rms_change, kurtosis_change, residual_rms_values, residual_kurtosis_values, mesh_harmonic_amp, mesh_sideband_amp, ...
    'VariableNames', {'CrackDepth_mm', 'RMS_m', 'Kurtosis', 'RMSChange_pct', 'KurtosisChange_pct', 'HealthySubtractedResidualRMS_m', 'HealthySubtractedResidualKurtosis', 'FourthMeshHarmonicAmp', 'FourthMeshLowerSidebandAmp'});
writetable(summary, fullfile(out_dir, 'summary_constant_depth.csv'));
end

function write_time_mesh_csv(out_dir, params, geom, crack_depths, time_mesh_profiles)
time = linspace(0, geom.t_rotation, geom.n_rotation)';
theta_deg = linspace(0, 360, geom.n_rotation)';
data = table(time, theta_deg, 'VariableNames', {'Time_s', 'PinionAngle_deg'});
for i = 1:numel(crack_depths)
    name = sprintf('K_q0_%03d_um_N_per_m', round(crack_depths(i) * 1e6));
    data.(name) = time_mesh_profiles(i, :)';
end
writetable(data, fullfile(out_dir, 'time_mesh_stiffness_rotational_period.csv'));
end

function make_figures(out_dir, params, geom, crack_depths, single_profiles, total_profiles, time_mesh_profiles, response_cases, rms_change, kurtosis_change, residual_rms_values, residual_kurtosis_values, mesh_harmonic_amp, mesh_sideband_amp)
phase_deg = linspace(0, 360 / params.N_p, geom.n_mesh);
rotation_time = linspace(0, geom.t_rotation, geom.n_rotation);
rotation_angle = linspace(0, 360, geom.n_rotation);
depth_labels = compose('q0 = %.1f mm', crack_depths * 1e3);

fig = figure('Visible', 'off');
plot(phase_deg, single_profiles' / 1e8, 'LineWidth', 1.2);
grid on
xlabel('Pinion angular position over one tooth pitch (deg)')
ylabel('Single-pair mesh stiffness (10^8 N/m)')
title('Single-pair stiffness, constant through-width crack')
legend(depth_labels, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'single_pair_stiffness.png'));
close(fig)

fig = figure('Visible', 'off');
plot(phase_deg, total_profiles' / 1e8, 'LineWidth', 1.2);
grid on
xlabel('Pinion angular position over one tooth pitch (deg)')
ylabel('Total mesh stiffness (10^8 N/m)')
title('Total time-varying mesh stiffness')
legend(depth_labels, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'total_mesh_stiffness.png'));
close(fig)

fig = figure('Visible', 'off');
plot(rotation_time, time_mesh_profiles' / 1e8, 'LineWidth', 1.0);
grid on
xlabel('Time over one pinion revolution (s)')
ylabel('Time-varying mesh stiffness (10^8 N/m)')
title('Mesh stiffness over one rotational period')
legend(depth_labels, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'time_mesh_stiffness_rotational_period.png'));
close(fig)

fig = figure('Visible', 'off');
plot(rotation_angle, time_mesh_profiles' / 1e8, 'LineWidth', 1.0);
grid on
xlabel('Pinion angle over one revolution (deg)')
ylabel('Time-varying mesh stiffness (10^8 N/m)')
title('Mesh stiffness over one pinion revolution')
legend(depth_labels, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'angle_mesh_stiffness_rotational_period.png'));
close(fig)

fig = figure('Visible', 'off');
hold on
for i = 1:numel(response_cases)
    plot(response_cases(i).t, response_cases(i).y_p * 1e6, 'LineWidth', 1.0);
end
hold off
grid on
xlabel('Time (s)')
ylabel('Pinion y displacement (um)')
title('Steady-state pinion displacement')
xlim([0.22, 0.30])
legend(depth_labels, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'pinion_y_displacement.png'));
close(fig)

fig = figure('Visible', 'off');
plot(crack_depths * 1e3, rms_change, 'o-', 'LineWidth', 1.2);
hold on
plot(crack_depths * 1e3, kurtosis_change, 's-', 'LineWidth', 1.2);
hold off
grid on
xlabel('Crack depth q0 (mm)')
ylabel('Change from healthy (%)')
title('RMS and kurtosis change')
legend({'RMS', 'Kurtosis'}, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'rms_kurtosis_change.png'));
close(fig)

fig = figure('Visible', 'off');
yyaxis left
plot(crack_depths * 1e3, residual_rms_values, 'o-', 'LineWidth', 1.2);
ylabel('Residual RMS (m)')
yyaxis right
plot(crack_depths * 1e3, residual_kurtosis_values, 's-', 'LineWidth', 1.2);
ylabel('Residual kurtosis')
grid on
xlabel('Crack depth q0 (mm)')
title('Healthy-subtracted residual indicators')
legend({'Residual RMS', 'Residual kurtosis'}, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'residual_indicators.png'));
close(fig)

fig = figure('Visible', 'off');
hold on
for i = 1:numel(response_cases)
    plot(response_cases(i).freq, response_cases(i).amp, 'LineWidth', 1.0);
end
hold off
grid on
xlabel('Frequency (Hz)')
ylabel('Amplitude (m)')
title('Pinion y displacement spectrum')
xlim([0, 10000])
legend(depth_labels, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'frequency_spectrum.png'));
close(fig)

fig = figure('Visible', 'off');
plot(crack_depths * 1e3, mesh_harmonic_amp, 'o-', 'LineWidth', 1.2);
hold on
plot(crack_depths * 1e3, mesh_sideband_amp, 's-', 'LineWidth', 1.2);
hold off
grid on
xlabel('Crack depth q0 (mm)')
ylabel('Amplitude (m)')
title('Fourth mesh harmonic and lower sideband')
legend({'4x mesh harmonic', '4x mesh lower sideband'}, 'Location', 'best')
saveas(fig, fullfile(out_dir, 'sideband_marker_change.png'));
close(fig)
end

function root = bisect_scalar(fun, a, b, tol)
fa = fun(a);
fb = fun(b);
if fa == 0
    root = a;
    return
end
if fb == 0
    root = b;
    return
end
if fa * fb > 0
    error('Bisection interval does not bracket a root: f(a)=%g, f(b)=%g', fa, fb);
end

root = (a + b) / 2;
fr = fun(root);
while abs(fr) > tol && abs(b - a) > tol
    if fa * fr < 0
        b = root;
        fb = fr;
    else
        a = root;
        fa = fr;
    end
    root = (a + b) / 2;
    fr = fun(root);
end
end

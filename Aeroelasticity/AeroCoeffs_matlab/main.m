clear; close all; clc;

addpath('functions/');
%% inizializzazione
[data, theta_target] = initialize_rotor_model('rigid'); 
N = data.N;

%% scalatura
K_blade_block = data.K_nom(data.idx_blades, data.idx_blades);
k_ref = mean(abs(diag(K_blade_block)));

scale.d0 = 1e-3 * k_ref;      
scale.ar = 1e-3 * k_ref;         
scale.ai = 1e-3 * k_ref;          
data.J_normalization = 1 / sum(sum(abs(data.y).^2));
scale.f  = 0.1;      
data.scale = scale;

%% guess primo tentativo
perturbation_magnitude = 0.05; 
d0_exact_norm = theta_target(1:N) / scale.d0;

% Initial Guess Completa
theta0_d0 = ones(N, 1);
theta0_ar = ones(N, 1);
theta0_ai = ones(N, 1); 

% Forcing Guess (EO=6)
Target_EO = 6; 
phase_vec = (0:N-1)' * (2*pi * Target_EO / N);
theta0_fr = 1.0 * cos(phase_vec);
theta0_fi = 1.0 * sin(phase_vec);

theta0_full = [theta0_d0; theta0_ar; theta0_ai; theta0_fr; theta0_fi];

% Bounds Completi
lb_d0 = -10 * ones(N, 1);    ub_d0 =  10 * ones(N, 1);
lb_ar = -10 * ones(N, 1);    ub_ar =  10 * ones(N, 1);
lb_ai = -20 * ones(N, 1);    ub_ai =  20 * ones(N, 1);
lb_fr = -20 * ones(N, 1);    ub_fr =  20 * ones(N, 1); 
lb_fi = -10 * ones(N, 1);    ub_fi =  10 * ones(N, 1);

lb_full = [lb_d0; lb_ar; lb_ai; lb_fr; lb_fi]; 
ub_full = [ub_d0; ub_ar; ub_ai; ub_fr; ub_fi]; 

%% Vincolo ar
% [d0 (1:N), ar (N+1:2N), ai, fr, fi]
idx_fixed = 2*N; %
fixed_val = 0;

idx_mask = true(5*N, 1);
idx_mask(idx_fixed) = false;

theta0_opt = theta0_full(idx_mask);
lb_opt     = lb_full(idx_mask);
ub_opt     = ub_full(idx_mask);

wrapper_fun = @(t) rotor_objective_wrapper_fixed(t, idx_mask, idx_fixed, fixed_val, data);

%% solver
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ... 
    'HessianApproximation', 'lbfgs', ...
    'Display', 'iter-detailed', ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 100000, ...
    'StepTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);

[theta_opt_reduced, final_J, exitflag] = fmincon(wrapper_fun, ...
                                    theta0_opt, [], [], [], [], lb_opt, ub_opt, [], options);

%% plot
theta_full_reconstructed = zeros(length(idx_mask), 1);
theta_full_reconstructed(idx_mask) = theta_opt_reduced;
theta_full_reconstructed(idx_fixed) = fixed_val;

d0_opt      = theta_full_reconstructed(1:N);
ar_opt      = theta_full_reconstructed(N+1 : 2*N);
ai_opt      = theta_full_reconstructed(2*N+1 : 3*N);
fr_opt      = theta_full_reconstructed(3*N+1 : 4*N);
fi_opt      = theta_full_reconstructed(4*N+1 : 5*N);

% variabili Target normalizzate
d0_tgt      = theta_target(1:N)           / scale.d0;
ar_tgt      = theta_target(N+1 : 2*N)     / scale.ar;
ai_tgt      = theta_target(2*N+1 : 3*N)   / scale.ai;
fr_tgt      = theta_target(3*N+1 : 4*N)   / scale.f;
fi_tgt      = theta_target(4*N+1 : 5*N)   / scale.f;

% Calcolo FFT Forzanti
f_vec_tgt_c = fr_tgt + 1i * fi_tgt;
f_vec_opt_c = fr_opt + 1i * fi_opt;
f_vec_guess_c = theta0_fr + 1i * theta0_fi;

F_nodal_tgt = fft(f_vec_tgt_c)/N;
F_nodal_opt = fft(f_vec_opt_c)/N;
F_nodal_guess = fft(f_vec_guess_c)/N;

sigma_val = data.sigmaNoise; 
err_d0 = ones(size(d0_opt)) * sigma_val;
err_ar = ones(size(ar_opt)) * sigma_val;
err_ai = ones(size(ai_opt)) * sigma_val;
err_F  = ones(N, 1) * (sigma_val * 0.5); 

if mod(N, 2) == 0, nd_axis = (-N/2 : N/2-1); else, nd_axis = (-(N-1)/2 : (N-1)/2); end

figure('Name', 'Results', 'Color', 'w', 'WindowStyle', 'docked');

col_bar = [1 1 0.2]; 
col_err = 'r';      
col_guess = 'g';    

% --- 1. Mistuning (d0) ---
ax1 = subplot(2, 2, 1);
hold on;

h_exact = bar(1:N, d0_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]);
h_guess = plot(1:N, theta0_d0, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', col_guess, 'Color', col_guess);
h_est   = errorbar(1:N, d0_opt, err_d0, 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 4, ...
    'MarkerFaceColor', col_err, 'Color', col_err);

h_fix_leg = plot(NaN, NaN, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 7, ...
                   'MarkerFaceColor', 'b', 'Color', 'b');
hold off;
grid on; title('Mistuning (d_0)'); 
xlabel('Blade Number'); ylabel('d_0'); xlim([0 N+1]);

% --- 2. Aero Stiffness (ar) ---
subplot(2, 2, 2);
bar(1:N, ar_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;

plot(1:N-1, theta0_ar(1:end-1), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', col_guess, 'Color', col_guess);

local_idx_fixed = idx_fixed - N; 

errorbar(1:N, ar_opt, err_ar, 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 4, ...
    'MarkerFaceColor', col_err, 'Color', col_err);
plot(local_idx_fixed, ar_opt(local_idx_fixed), 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'Color', 'b', 'LineWidth', 0.5);

hold off;
grid on; title('Aero Stiffness (Real)');
xlabel('Blade Number'); ylabel('a_02'); xlim([0 N+1]);

% --- 3. Aero Damping (ai) ---
subplot(2, 2, 3);
bar(1:N, ai_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;
plot(1:N, theta0_ai, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', col_guess, 'Color', col_guess);
errorbar(1:N, ai_opt, err_ai, 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 4, ...
    'MarkerFaceColor', col_err, 'Color', col_err);
hold off;
grid on; title('Aero Damping (Imag)');
xlabel('Blade Number'); ylabel('a_1'); xlim([0 N+1]);

% --- 4. Forcing Spectrum FFT ---
subplot(2, 2, 4);
bar(nd_axis, abs(F_nodal_tgt), 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;
plot(nd_axis, abs(F_nodal_guess), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', col_guess, 'Color', col_guess);
errorbar(nd_axis, abs(F_nodal_opt), err_F, 'LineStyle', 'none', 'Marker', 's', 'MarkerSize', 4, ...
    'MarkerFaceColor', col_err, 'Color', col_err);
hold off;
xlabel('Nodal Diameter'); title('Forcing Spectrum'); 
grid on; xlim([-16 16]);

% --- LEGENDA ESTERNA MANUALE ---

lgd = legend([h_exact, h_guess, h_est, h_fix_leg], ...
       {'Exact Model', 'Initial Guess', 'ML Estimate', 'Fixed Param'}, ...
       'Orientation', 'horizontal');
set(lgd, 'Position', [0.20, 0.94, 0.7, 0.05]); 
set(lgd, 'Box', 'off', 'Color', 'w');


% --- Figura 2: Response Graph (DOCKED) ---
d0_phys_opt = d0_opt * scale.d0;
ar_phys_opt = ar_opt * scale.ar;
ai_phys_opt = ai_opt * scale.ai;
F_phys_opt  = (fr_opt + 1i * fi_opt) * scale.f;

d0_phys_tgt = d0_tgt * scale.d0;
ar_phys_tgt = ar_tgt * scale.ar;
ai_phys_tgt = ai_tgt * scale.ai;
F_phys_tgt  = (fr_tgt + 1i * fi_tgt) * scale.f;

if isfield(data, 'omega')
    w_vec = data.omega;
else
    w_vec = linspace(16000, 18000, 500); 
end
K_nom = data.K_nom;
M = data.M_mat;

n_freq = length(w_vec);
amp_model = zeros(1, n_freq);
amp_ident = zeros(1, n_freq);
amp_tuned = zeros(1, n_freq);

for k = 1:n_freq
    w = w_vec(k);
    w2 = w^2;
    Z_base = K_nom - w2 * M;
    
    % Model
    Z_model = Z_base + diag(d0_phys_tgt) + diag(ar_phys_tgt) + 1i * diag(ai_phys_tgt);
    amp_model(k) = max(abs(Z_model \ F_phys_tgt));
    
    % Identified
    Z_ident = Z_base + diag(d0_phys_opt) + diag(ar_phys_opt) + 1i * diag(ai_phys_opt);
    amp_ident(k) = max(abs(Z_ident \ F_phys_opt));
    
    % Tuned (Z_tuned=Z_base, F_tuned=F_phys_opt)
    Z_tuned = Z_base; 
    amp_tuned(k) = max(abs(Z_tuned \ F_phys_opt));
end

figure('Name', 'Log Frequency Response', 'Color', 'w', 'WindowStyle', 'docked');
hold on;
amp_model_log = log10(abs(amp_model));
amp_ident_log = log10(abs(amp_ident));
amp_tuned_log = log10(abs(amp_tuned));
w_vec_hz = w_vec / (2*pi);

if isfield(data, 'y') && ~isempty(data.y)
    plot(w_vec_hz, log10(abs(data.y)), 'Color', [0.85 0.85 0.85], ...
        'LineWidth', 0.5, 'HandleVisibility', 'off');
end

h_mod = plot(w_vec_hz, amp_model_log, 'k-', 'LineWidth', 2, 'DisplayName', 'Model System');
h_ide = plot(w_vec_hz, amp_ident_log, 'r--', 'LineWidth', 2, 'DisplayName', 'Identified System');
h_tun = plot(w_vec_hz, amp_tuned_log, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Tuned System');

grid on; grid minor; 
xlabel('Frequency [Hz]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Vibration Amplitude [m] - Log Scale', 'FontSize', 11, 'FontWeight', 'bold');
title('Forced Response Reconstruction (Log Scale)', 'FontSize', 12);
xlim([min(w_vec_hz) max(w_vec_hz)]);

valid_vals = [amp_model_log, amp_ident_log, amp_tuned_log];
valid_vals = valid_vals(~isinf(valid_vals));
if ~isempty(valid_vals)
    max_val = max(valid_vals);
    min_val = min(valid_vals);
    ylim([min_val, max_val + 0.5]);
end
legend([h_mod, h_ide, h_tun], 'Location', 'northeast');
box on;
hold off;

%% Helper function 
function [J, grad_reduced] = rotor_objective_wrapper_fixed(theta_reduced, idx_mask, idx_fixed, fixed_val, data)
    theta_full = zeros(length(idx_mask), 1);
    theta_full(idx_mask) = theta_reduced;
    theta_full(idx_fixed) = fixed_val;
    
    [J, grad_full] = rotor_objective_with_grad(theta_full, data);
    
    grad_reduced = grad_full(idx_mask);
end

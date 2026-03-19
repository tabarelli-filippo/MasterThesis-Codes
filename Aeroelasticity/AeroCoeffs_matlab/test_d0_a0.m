clear; close all; clc;

%% Inizializzazione Ambiente e Dati
[data, theta_target] = initialize_rotor_model('rigid'); 
N = data.N;

%% definizione fattori scala
K_blade_block = data.K_nom(data.idx_blades, data.idx_blades);
k_ref = mean(abs(diag(K_blade_block)));

scale.d0 = 1e-3*k_ref ;      % Mistuning: k_phys = d0_norm * k_ref
scale.ar = 1e-5*k_ref ;          % Aero Real:  Re(A)_phys = ar_norm * k_ref
scale.ai = 1e-4*k_ref ;          % Aero Imag:  Im(A)_phys = ai_norm * k_ref
data.J_normalization = 1 / sum(sum(abs(data.y).^2));
scale.f  = 1.0;      

data.scale = scale;

%% bounds
tgt_d0_norm = theta_target(1 : N) / scale.d0;
tgt_ar_norm = theta_target(N+1 : 2*N) / scale.ar;
tgt_ai_norm = theta_target(2*N+1 : 3*N) / scale.ai;
tgt_fr_norm = theta_target(3*N+1 : 4*N) / scale.f;
tgt_fi_norm = theta_target(4*N+1 : 5*N) / scale.f;

theta0_d0 = ones(N,1);
theta0_ar = ones(N-1, 1);  

lb_d0 = -10 * ones(N, 1);   ub_d0 =  10 * ones(N, 1);
lb_ar = -10 * ones(N-1, 1);   ub_ar =  10 * ones(N-1, 1);

lb = [lb_d0; lb_ar];       ub = [ub_d0; ub_ar];

theta0 = [theta0_d0; theta0_ar];

%  Wrapper function for d0 and ar
fixed_params = [tgt_ai_norm; tgt_fr_norm; tgt_fi_norm];
fixed_ar_val = 0;

% costraint on d0 - sum of d0 = 0
Aeq = [];   
beq = [];

wrapper_obj = @(t) rotor_objective_wrapper(t, fixed_params, fixed_ar_val, data);

%% Global minimum
% rng default 
% num_start_points = 2*N; 
% 
% fprintf('Generazione di %d punti tramite Ipercubo Latino (Maximin)...\n', num_start_points);
% lhs_unit = lhsdesign(num_start_points, 2*N-1, 'criterion', 'maximin', 'iterations', 50);
% 
% start_points_matrix = repmat(lb', num_start_points, 1) + ...
%                       repmat((ub - lb)', num_start_points, 1) .* lhs_unit;
% 
% tpoints = CustomStartPointSet(start_points_matrix);
% 
% local_opts = optimoptions(@fmincon, ...
%     'Algorithm', 'interior-point', ...
%     'SpecifyObjectiveGradient', true, ... 
%     'HessianApproximation', 'bfgs', ...
%     'OptimalityTolerance', 1e-6, ...
%     'StepTolerance', 1e-8, ...
%     'MaxIterations', 1000, ...
%     'Display', 'off'); 
% 
% problem = createOptimProblem('fmincon', ...
%     'objective', wrapper_obj, ...
%     'x0', theta0, ...
%     'lb', lb, ...
%     'ub', ub, ...
%     'options', local_opts);
% 
% ms = MultiStart(...
%     'Display', 'iter', ...
%     'UseParallel', true); 
% 
% fprintf('Avvio MultiStart con %d punti di partenza...\n', num_start_points);
% [startPoint, f, exitflag1, output1, solutions] = run(ms, problem, tpoints);
% 

startPoint = theta0;
%% solver
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ...
    'HessianApproximation', 'lbfgs', ...
    'Display', 'iter-detailed', ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 20000, ...
    'StepTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);

[theta_norm_opt, final_J] = fmincon(wrapper_obj, ...
                                    startPoint, [], [], Aeq, beq, lb, ub, [], options);

%% grafici
theta_full_reconstructed = [theta_norm_opt(1:2*N-1); fixed_ar_val; fixed_params];

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

sigma_val = data.sigmaNoise; 
err_d0 = ones(size(d0_opt)) * sigma_val;
err_ar = ones(size(ar_opt)) * sigma_val;
err_ai = ones(size(ai_opt)) * sigma_val;
err_F  = ones(N, 1) * (sigma_val * 0.5); 

figure('Name', 'Results', 'Color', 'w');
% Colors
col_bar = [1 1 0.2]; % Yellow (Exact Model)
col_err = 'r';       % Red (Identified)

% 1. Mistuning
subplot(2, 2, 1);
bar(1:N, d0_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5], 'DisplayName', 'Exact Model'); hold on;
hold on;

plot(1:N, theta0_d0, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'g', ...
    'Color', 'g', 'LineWidth', 1, ...
    'DisplayName', 'Initial Guess');

errorbar(1:N, d0_opt, err_d0, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1, ...
    'DisplayName', 'ML Estimate');
hold off;

grid on; legend('Location','best'); title('Mistuning (d_0)'); 
xlabel('Blade Number');
ylabel('d_0')
xlim([0 N+1]);

% 2. Aero Stiffness
subplot(2, 2, 2);
bar(1:N, ar_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;
hold on;
plot(1:N-1, theta0_ar, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'g', ...
    'Color', 'g', 'LineWidth', 1, ...
    'DisplayName', 'Initial Guess');
errorbar(1:N, ar_opt, err_ar, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1, ...
    'DisplayName', 'ML Estimate');

hold off;
grid on; title('Aero Stiffness (Real)');
xlabel('Blade Number');
ylabel('a_02')
xlim([0 N+1]);

% 3. Aero Damping
subplot(2, 2, 3);
bar(1:N, ai_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;
errorbar(1:N, ai_opt, err_ai, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1);
grid on; title('Aero Damping (Imag)');
xlabel('Blade Number');
ylabel('a_1')
xlim([0 N+1]);

% 4. Forcing Spectrum FFT
f_vec_tgt_c = fr_tgt + 1i * fi_tgt;
f_vec_opt_c = fr_opt + 1i * fi_opt;

F_nodal_tgt = fft(f_vec_tgt_c)/N;
F_nodal_opt = fft(f_vec_opt_c)/N;

if mod(N, 2) == 0, nd_axis = (-N/2 : N/2-1); else, nd_axis = (-(N-1)/2 : (N-1)/2); end

subplot(2, 2, 4);
bar(nd_axis, abs(F_nodal_tgt), 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;
errorbar(nd_axis, abs(F_nodal_opt), err_F, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1);
xlabel('Nodal Diameter'); title('Forcing Spectrum'); 
grid on; xlim([-16 16]);

% Response Graph
d0_phys_opt = d0_opt * scale.d0;
ar_phys_opt = ar_opt * scale.ar;
ai_phys_opt = ai_opt * scale.ai;
F_phys_opt  = (fr_opt + 1i * fi_opt) * scale.f;

d0_phys_tgt = d0_tgt * scale.d0;
ar_phys_tgt = ar_tgt * scale.ar;
ai_phys_tgt = ai_tgt * scale.ai;
F_phys_tgt  = (fr_tgt + 1i * fi_tgt) * scale.f;

d0_phys_tuned = zeros(N, 1); 
ar_phys_tuned = ar_phys_opt; 
ai_phys_tuned = ai_phys_opt;
F_phys_tuned  = F_phys_opt;

if isfield(data, 'omega')
    w_vec = data.omega;
else
    w_vec = linspace(16000, 18000, 500); 
end

K_nom = data.K_nom;
M = data.M_mat;

% Pre-allocation
n_freq = length(w_vec);
amp_model = zeros(1, n_freq);
amp_ident = zeros(1, n_freq);
amp_tuned = zeros(1, n_freq);

for k = 1:n_freq
    w = w_vec(k);
    w2 = w^2;
    
    Z_base = K_nom - w2 * M;
    
    % Model System
    Z_model = Z_base + diag(d0_phys_tgt) + diag(ar_phys_tgt) + 1i * diag(ai_phys_tgt);
    X_model = Z_model \ F_phys_tgt;
    amp_model(k) = max(abs(X_model));
    
    % Identified System
    Z_ident = Z_base + diag(d0_phys_opt) + diag(ar_phys_opt) + 1i * diag(ai_phys_opt);
    X_ident = Z_ident \ F_phys_opt;
    amp_ident(k) = max(abs(X_ident));
    
    % Tuned System
    Z_tuned = Z_base;
    X_tuned = Z_tuned \ F_phys_tuned;
    amp_tuned(k) = max(abs(X_tuned));
end

figure(); clf;
set(gcf, 'Name', 'Logarithmic Frequency Response ', 'Color', 'w', 'WindowStyle', 'docked');
hold on;

amp_model(k) = max(abs(X_model));
amp_ident(k) = max(abs(X_ident));
amp_tuned(k) = max(abs(X_tuned));

amp_model = log10(abs(amp_model));
amp_ident = log10(abs(amp_ident));
amp_tuned = log10(abs(amp_tuned));

w_vec = w_vec /2/pi;

if isfield(data, 'y') && ~isempty(data.y)
    h_raw = plot(w_vec, log10(abs(data.y)), 'Color', [0.85 0.85 0.85], ...
        'LineWidth', 0.5, 'HandleVisibility', 'off');
end

h_mod = plot(w_vec, amp_model, 'k-', 'LineWidth', 2, 'DisplayName', 'Model System');
h_ide = plot(w_vec, amp_ident, 'r--', 'LineWidth', 2, 'DisplayName', 'Identified System');
h_tun = plot(w_vec, amp_tuned, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Tuned System');

grid on; grid minor; 

xlabel('Frequency [Hz]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Vibration Amplitude [m] - Log Scale', 'FontSize', 11, 'FontWeight', 'bold');
title('Forced Response Reconstruction (Log Scale)', 'FontSize', 12);

xlim([min(w_vec) max(w_vec)]);

max_val = max([max(amp_model), max(amp_ident), max(amp_tuned)]);
min_val = min([min(amp_model), min(amp_ident), min(amp_tuned)]);
ylim([min_val, max_val*0.8]);

legend([h_mod, h_ide, h_tun], 'Location', 'northeast');
box on;
hold off;

%% Helper function
function [J, grad] = rotor_objective_wrapper(vars, fixed_params,fixed_ar_val , data)
    theta_full = [vars; fixed_ar_val; fixed_params];
    
    [J, grad_full] = rotor_objective_with_grad(theta_full, data);
    N = data.N;
    grad = grad_full(1:2*N-1);
end

clear; close all; clc;

addpath('functions/');
%% Inizializzazione Ambiente e Dati
[data, theta_target] = initialize_rotor_model('rigid'); 
N = data.N;

%% definizione fattori scala
K_blade_block = data.K_nom(data.idx_blades, data.idx_blades);
k_ref = mean(abs(diag(K_blade_block)));

scale.d0 = 1e-3* k_ref ;      % Mistuning: k_phys = d0_norm * k_ref
scale.ar = 1e-5* k_ref ;          % Aero Real:  Re(A)_phys = ar_norm * k_ref
scale.ai = 1e-4* k_ref ;          % Aero Imag:  Im(A)_phys = ai_norm * k_ref
data.J_normalization = 1 / sum(sum(abs(data.y).^2));
scale.f  = 0.1;      

data.scale = scale;

%% bounds
tgt_d0_norm = theta_target(1 : N) / scale.d0;
tgt_ar_norm = theta_target(N+1 : 2*N) / scale.ar;
tgt_fr_norm = theta_target(3*N+1 : 4*N) / scale.f;
tgt_fi_norm = theta_target(4*N+1 : 5*N) / scale.f;

theta0_ai = ones(N, 1); 

lb_ai = -10 * ones(N, 1);
ub_ai =  10 * ones(N, 1);

%  Wrapper function for mistuning only variables
fixed_params1 = [tgt_d0_norm; tgt_ar_norm];
fixed_params2 = [tgt_fr_norm; tgt_fi_norm];
wrapper_obj = @(t) rotor_objective_wrapper(t, fixed_params1, fixed_params2, data);

%% solver
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient', true, ... 
    'HessianApproximation', 'lbfgs', ...
    'Display', 'iter-detailed', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 10000, ...
    'StepTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-8);

[theta_norm_opt, final_J] = fmincon(wrapper_obj, ...
                                    theta0_ai, [], [], [], [], lb_ai, ub_ai, [], options);

%% grafici
theta_full_reconstructed = [fixed_params1; theta_norm_opt; fixed_params2];

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

figure('Name', 'Results', 'WindowStyle', 'docked', 'Color', 'w');

% Colors
col_bar = [1 1 0.2]; % Yellow (Exact Model)
col_err = 'r';       % Red (Identified)

% 1. Mistuning
subplot(2, 2, 1);
bar(1:N, d0_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5], 'DisplayName', 'Exact Model'); hold on;
errorbar(1:N, d0_opt, err_d0, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1, 'DisplayName', 'ML Estimate');
grid on; legend('Location','best'); title('Mistuning (d_0)'); 
xlabel('Blade Number');
ylabel('d_0')
xlim([0 N+1]);

% 2. Aero Stiffness
subplot(2, 2, 2);
bar(1:N, ar_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;
errorbar(1:N, ar_opt, err_ar, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1);
grid on; title('Aero Stiffness (Real)');
xlabel('Blade Number');
ylabel('a_0')
xlim([0 N+1]);

% 3. Aero Damping
subplot(2, 2, 3);
bar(1:N, ai_tgt, 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5], 'DisplayName', 'Exact Model'); 
hold on;
plot(1:N, theta0_ai, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'g', ...
    'Color', 'g', 'LineWidth', 1, ...
    'DisplayName', 'Initial Guess');
errorbar(1:N, ai_opt, err_ai, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1, ...
    'DisplayName', 'ML Estimate');
hold off;
legend('Location','best');
grid on; title('Aero Damping (Imag)');
xlabel('Blade Number');
ylabel('a_1')
xlim([0 N+1]);

% 4. Forcing Spectrum FFT
f_vec_tgt_c = fr_tgt + 1i * fi_tgt;
f_vec_opt_c = fr_opt + 1i * fi_opt;

F_nodal_tgt = fftshift(fft(f_vec_tgt_c)/N);
F_nodal_opt = fftshift(fft(f_vec_opt_c)/N);

if mod(N, 2) == 0, nd_axis = (-N/2 : N/2-1); else, nd_axis = (-(N-1)/2 : (N-1)/2); end

subplot(2, 2, 4);
bar(nd_axis, abs(F_nodal_tgt), 'FaceColor', col_bar, 'EdgeColor', [0.5 0.5 0.5]); hold on;
errorbar(nd_axis, abs(F_nodal_opt), err_F, 'LineStyle', 'none', ...
    'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', col_err, ...
    'Color', col_err, 'LineWidth', 1);
xlabel('Nodal Diameter'); title('Forcing Spectrum'); 
grid on; xlim([-16 16]);


function [J, grad_ai] = rotor_objective_wrapper(theta, fixed_params1, fixed_params2, data)
    theta_full = [fixed_params1; theta; fixed_params2];
    
    [J, grad_full] = rotor_objective_with_grad(theta_full, data);
    N = data.N;
    grad_ai = grad_full(2*N+1 : 3*N);
end

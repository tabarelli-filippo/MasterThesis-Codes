clear; close all; clc;

%% 1. Inizializzazione Ambiente e Dati
[data, theta_target] = initialize_rotor_model('rigid'); 
N = data.N;

%% 2. Definizione Fattori di Scala (NORMALIZZAZIONE)
K_blade_block = data.K_nom(data.idx_blades, data.idx_blades);
k_ref = mean(diag(K_blade_block));

scale.d0 = 0.05*k_ref ;      % Mistuning: k_phys = d0_norm * k_ref
scale.ar = k_ref ;      % Aero Real:  Re(A)_phys = ar_norm * k_ref
scale.ai = k_ref ;      % Aero Imag:  Im(A)_phys = ai_norm * k_ref
data.J_normalization = 1 / sum(sum(abs(data.y).^2));
scale.f  = 10.0;      

data.scale = scale;

%% 3. Definizione Bounds
tgt_ar_norm = theta_target(N+1 : 2*N) / scale.ar;
tgt_ai_norm = theta_target(2*N+1 : 3*N) / scale.ai;
tgt_fr_norm = theta_target(3*N+1 : 4*N) / scale.f;
tgt_fi_norm = theta_target(4*N+1 : 5*N) / scale.f;

% Bounds per Mistuning (Lasciamo libertà al solver)
lb_d0 = -5 * ones(N, 1); 
ub_d0 =  5 * ones(N, 1);

% Bounds per parametri fissi (Usa un epsilon piccolissimo per stabilità numerica)
eps_b = 1e-5;
lb_ar = tgt_ar_norm - eps_b; ub_ar = tgt_ar_norm + eps_b;
lb_ai = tgt_ai_norm - eps_b; ub_ai = tgt_ai_norm + eps_b;
lb_fr  = tgt_fr_norm - eps_b; ub_fr  = tgt_fr_norm + eps_b;
lb_fi  = tgt_fi_norm - eps_b; ub_fi  = tgt_fi_norm + eps_b;

lb = [lb_d0; lb_ar; lb_ai; lb_fr; lb_fi];
ub = [ub_d0; ub_ar; ub_ai; ub_fr; ub_fi];

% Initial Guess 
theta0 = zeros(size(lb));
theta0(1:N) = 0; % Partiamo da zero mistuning
theta0(N+1 : 2*N)   = tgt_ar_norm;
theta0(2*N+1 : 3*N) = tgt_ai_norm;
theta0(3*N+1 : 4*N) = tgt_fr_norm;
theta0(4*N+1 : 5*N) = tgt_fi_norm;

figure('Name', 'Debug: Theta0 vs Bounds', 'Color', 'w', 'WindowStyle', 'docked');

% Definizioni indici per slicing del vettore
idx_d0 = 1:N;
idx_ar = N+1 : 2*N;
idx_ai = 2*N+1 : 3*N;
idx_f  = 3*N+1 : 5*N; % Include sia parte Reale che Immaginaria

% Helper function per plot rapido
plot_bounds_check(2, 2, 1, idx_d0, theta0, lb, ub, 'Mistuning (d_0)');
plot_bounds_check(2, 2, 2, idx_ar, theta0, lb, ub, 'Aero Real (a_r)');
plot_bounds_check(2, 2, 3, idx_ai, theta0, lb, ub, 'Aero Imag (a_i)');
plot_bounds_check(2, 2, 4, idx_f,  theta0, lb, ub, 'Forcing (Re & Im)');

%% Funzione Locale di Plotting
function plot_bounds_check(r, c, p, idx, t0, low, up, title_str)
    subplot(r, c, p); hold on; grid on; box on;
    
    % Estrazione dati locali
    val = t0(idx);
    l   = low(idx);
    u   = up(idx);
    x   = 1:length(idx);
    
    % 1. Plot Area Valida (Opzionale, visivamente utile)
    % Crea un poligono grigio chiaro tra lb e ub
    fill([x, fliplr(x)], [l', fliplr(u')], [0.9 0.9 0.9], 'EdgeColor', 'none');
    
    % 2. Plot Bounds
    plot(x, l, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Lower Bound');
    plot(x, u, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Upper Bound');
    
    % 3. Identifica Violazioni e Clamping
    tol = 1e-6;
    is_violation = (val < l - tol) | (val > u + tol);
    is_clamped   = (abs(val - l) <= tol) | (abs(val - u) <= tol);
    is_safe      = ~(is_violation | is_clamped);
    
    % 4. Plot Valori Theta0
    % Punti sicuri (Verdi o Blu)
    if any(is_safe)
        plot(x(is_safe), val(is_safe), 'b.', 'MarkerSize', 12, 'DisplayName', '\theta_0 Safe');
    end
    
    % Punti Clamped (Gialli/Arancio - Attenzione: Solver bloccato sul bordo)
    if any(is_clamped)
        plot(x(is_clamped), val(is_clamped), 'o', 'MarkerSize', 6, ...
            'MarkerEdgeColor', [0.85 0.33 0.1], 'LineWidth', 1.5, 'DisplayName', 'Clamped');
    end
    
    % Punti Violazione (Rossi/Croci - Errore Fatale per Interior-Point)
    if any(is_violation)
        plot(x(is_violation), val(is_violation), 'rx', 'MarkerSize', 10, ...
            'LineWidth', 2, 'DisplayName', 'VIOLATION');
    end
    
    % Formatting
    title(title_str, 'FontWeight', 'bold');
    xlabel('Index'); ylabel('Normalized Value');
    xlim([0 length(idx)+1]);
    
    % Aggiusta Y-Limits per vedere bene i marker anche se fuori bound
    y_all = [l; u; val];
    margin = (max(y_all) - min(y_all)) * 0.1;
    if margin == 0, margin = 1; end
    ylim([min(y_all)-margin, max(y_all)+margin]);
    
    % Legend (solo primo plot per pulizia o smart location)
    if p == 1
        legend('Location', 'best');
    end
end
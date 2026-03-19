function [data, theta_true] = initialize_rotor_model(model_type)
%INITIALIZE_ROTOR_MODEL Genera dati calibrati sul Purdue 3-Stage Compressor (P3S)
%
% Output:
%   data        : Struttura dati contenente matrici, frequenze e misure
%   theta_true  : Vettore parametri veri [d0; ar; ai; fr; fi]

    if nargin < 1, model_type = 'rigid'; end
    fprintf('Inizializzazione Modello P3S: %s\n', upper(model_type));
    
    %% 1. Parametri Fisici del Rotore
    data.N = 38;            % Numero pale
    f0_Hz  = 2727.87;       % Frequenza naturale media (Modo 1T)
    w0     = 2*pi*f0_Hz;  
    
    % Sorgente di eccitazione
    engine_order = 33;      
    
    % Sweep di frequenza (Attorno alla risonanza)
    rpm_center = f0_Hz * 60 / engine_order; 
    rpm_span   = 200;       
    N_steps    = 1000; 
    
    rpm_vec = linspace(rpm_center - rpm_span, rpm_center + rpm_span, N_steps);
    data.omega = rpm_vec / 30 * pi * engine_order; % rad/s
    data.rpm   = rpm_vec;
    data.freq_Hz = data.omega / (2*pi);

    %% 2. Matrici del Sistema (K_nom, M_mat)
    rng(42); 
    
    % Definizione differenziata dello smorzamento strutturale
    zeta_disk  = 0.0001; % 0.01% Smorzamento strutturale solo sul Disco
    zeta_blade = 1e-6;   % Smorzamento strutturale Pale (trascurabile)
    
    factor_disk  = (1 + 1i * 2 * zeta_disk);
    factor_blade = (1 + 1i * 2 * zeta_blade);
    
    switch lower(model_type)
        case 'rigid'
            % Modello Lumped Parameter (Solo pale, nessun disco)
            data.M_mat = eye(data.N);
            % Si applica solo il fattore blade (trascurabile)
            data.K_nom = (w0^2) * factor_blade * eye(data.N);
            data.idx_blades = 1:data.N;
            
        case 'flexible'
            % Modello Disco + Pale
            M_bb = eye(data.N);
            % Smorzamento pale trascurabile
            K_bb = (w0^2) * factor_blade * eye(data.N);
            
            N_disk = data.N;
            M_dd = eye(N_disk);
            nd = linspace(-1, 1, N_disk); 
            freqs_disk = w0 * (0.8 + 0.4 * abs(nd).^2); 
            
            % Applicazione smorzamento 0.01% SOLO al blocco disco
            K_dd = diag(freqs_disk.^2) * factor_disk;
            
            coup_factor = 0.05; 
            M_bd = coup_factor * (rand(data.N, N_disk) - 0.5); 
            % Il coupling stiffness rimane reale (o con smorzamento base)
            K_bd = coup_factor * w0^2 * (rand(data.N, N_disk) - 0.5);
            
            data.M_mat = [M_dd, M_bd'; M_bd, M_bb];
            data.K_nom = [K_dd, K_bd'; K_bd, K_bb];
            data.idx_blades = (N_disk + 1) : (N_disk + data.N);
            
        otherwise
            error('Tipo modello non supportato.');
    end
    
    %% 3. Definizione Parametri Veri
    rng(1); 
    
    % A. Mistuning (d0)
    sigma_freq_percent = 0.02; 
    sigma_mistuning_k  = 2 * sigma_freq_percent * w0^2; 
    
    d0_true = sigma_mistuning_k * randn(data.N, 1);
    
    % B. Aerodinamica (Coefficienti a_real, a_imag)
    zeta_target = 1e-3; 
    k_damping_base = w0^2 * 1e-4; 
    
    k_coupling_base = w0^2 * 1e-5;
    a_real_true = zeros(data.N, 1);
    a_imag_true = zeros(data.N, 1);
    
    % Pattern Aerodinamico
    a_real_true(end)     =  0;                  
    a_real_true(end-1)     =  4 * k_coupling_base; 
    a_real_true(1)   =  4 * k_coupling_base; 
    
    a_imag_true(1)     =  1 * k_damping_base; 
    a_imag_true(end-1)     =  3 * k_damping_base; 
    a_imag_true(end)   =  0 * k_damping_base; 
    
    % C. Forzante (fr, fi)
    f_scalar_amp = 0.148; 
    
    % Onda viaggiante
    wave_shape = exp(1i * (0:data.N-1)' * (2*pi * engine_order / data.N));
    f_vec_complex = f_scalar_amp * wave_shape;
    
    f_real_true = real(f_vec_complex);
    f_imag_true = imag(f_vec_complex);
    
    % Assemblaggio theta_true
    theta_true = [d0_true; a_real_true; a_imag_true; f_real_true; f_imag_true];
    
    %% 4. Calcolo Risposta (Forward Run)
    mk_circ = @(v) toeplitz(v(end:-1:1), [v(end); v(1:end-1)]);
    
    A_aero_mat = mk_circ(a_real_true) + 1i * mk_circ(a_imag_true);
    D0_mat     = diag(d0_true);
    
    idx = data.idx_blades;
    
    % --- Sistema TUNED ---
    A_tuned = data.K_nom;
    A_tuned(idx, idx) = A_tuned(idx, idx) + A_aero_mat; 
    
    % --- Sistema MISTUNED ---
    A_mistuned = data.K_nom;
    A_mistuned(idx, idx) = A_mistuned(idx, idx) + D0_mat + A_aero_mat;
    
    % Vettore forza full
    F_sys = zeros(size(data.M_mat, 1), 1);
    F_sys(idx) = f_vec_complex;
    
    % Buffer risultati
    y_mistuned_complex = zeros(data.N, length(data.omega));
    y_tuned_complex    = zeros(data.N, length(data.omega));
    
    % Solver in frequenza
    for k = 1:length(data.omega)
        wk = data.omega(k);
        Dyn_Mat_Base = -(wk^2) * data.M_mat;
        
        % Risposta Tuned
        xt = (A_tuned + Dyn_Mat_Base) \ F_sys;
        y_tuned_complex(:, k) = xt(idx);
        
        % Risposta Mistuned
        xm = (A_mistuned + Dyn_Mat_Base) \ F_sys;
        y_mistuned_complex(:, k) = xm(idx);
    end
    
    %% 5. Output e Rumore
    noise_percent = 0.0; 
    max_amp = max(abs(y_mistuned_complex(:)));
    data.sigmaNoise = noise_percent * max_amp;
    
    noise_mat = data.sigmaNoise * (randn(size(y_mistuned_complex)) + 1i*randn(size(y_mistuned_complex)));
    data.y = y_mistuned_complex + noise_mat;
    
    data.y_tuned = y_tuned_complex; 
    data.y_clean = y_mistuned_complex; 
    
    data.scale.d0 = 1; data.scale.ar = 1; data.scale.ai = 1; data.scale.f = 1;
    data.J_normalization = 1;
    
    fprintf('Dati P3S generati. Freq: %.1f Hz, EO: %d, Smorzamento: %.2e, Forza: %.3f N\n', ...
        f0_Hz, engine_order, zeta_target, f_scalar_amp);
end
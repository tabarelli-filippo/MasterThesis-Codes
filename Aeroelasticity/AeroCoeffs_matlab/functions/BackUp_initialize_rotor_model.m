function [data, theta_true] = initialize_rotor_model(model_type)
%INITIALIZE_ROTOR_MODEL Genera dati compatibili con rotor_objective_with_grad e test_d0
%
% Output:
%   data        : Struttura dati contenente matrici, frequenze e misure (data.y complesso)
%   theta_true  : Vettore parametri veri [d0; ar; ai; fr; fi]

    if nargin < 1, model_type = 'rigid'; end
    fprintf('Inizializzazione Modello: %s (Compatibile con Solver e Plot)\n', upper(model_type));
    
    %% 1. Parametri Fisici del Rotore
    data.N = 33;          
    f0_Hz  = 2760;        
    w0     = 2*pi*f0_Hz;  
    
    engine_order = 6;     
    
    % Sweep di frequenza (Attorno alla risonanza)
    rpm_center =  f0_Hz * 60 / engine_order; 
    rpm_span   = 1200; 
    N_steps    = 500; % Numero punti per lo sweep
    
    rpm_vec = linspace(rpm_center - rpm_span, rpm_center + rpm_span, N_steps);
    data.omega = rpm_vec / 30 * pi * engine_order; % rad/s
    data.rpm   = rpm_vec;
    
    % Asse Frequenza in Hz (utile per i plot)
    data.freq_Hz = data.omega / (2*pi);

    %% 2. Matrici del Sistema (K_nom, M_mat)
    rng(42); % Seed per riproducibilità strutturale
    
    % Smorzamento strutturale (0.01% -> zeta = 0.0001)
    % Nota: rotor_objective_with_grad usa K_nom. 
    % Inseriamo lo smorzamento come parte immaginaria della rigidezza.
    zeta_struc = 0.0001; 
    complex_stiff_factor = (1 + 1i * 2 * zeta_struc);
    
    switch lower(model_type)
        case 'rigid'
            % Modello a 1 GDL per settore (solo pala)
            data.M_mat = eye(data.N);
            data.K_nom = (w0^2) * complex_stiff_factor * eye(data.N);
            data.idx_blades = 1:data.N;
            
        case 'flexible'
            % Modello Disco + Pale
            M_bb = eye(data.N);
            K_bb = (w0^2) * complex_stiff_factor * eye(data.N);
            
            N_disk = data.N;
            M_dd = eye(N_disk);
            nd = linspace(-1, 1, N_disk); 
            freqs_disk = w0 * (0.8 + 0.4 * abs(nd).^2); 
            
            K_dd = diag(freqs_disk.^2) * complex_stiff_factor;
            
            coup_factor = 0.05; 
            M_bd = coup_factor * (rand(data.N, N_disk) - 0.5); 
            K_bd = coup_factor * w0^2 * (rand(data.N, N_disk) - 0.5);
            
            data.M_mat = [M_dd, M_bd'; M_bd, M_bb];
            data.K_nom = [K_dd, K_bd'; K_bd, K_bb];
            % Indici delle pale (dopo i GDL del disco)
            data.idx_blades = (N_disk + 1) : (N_disk + data.N);
            
        otherwise
            error('Tipo modello non supportato.');
    end
    
    %% 3. Definizione Parametri Veri (Ground Truth)
    rng(1); % Seed per parametri fisici
    
    % A. Mistuning (d0)
    sigma_mistuning = 0.004 * w0^2; 
    d0_true = sigma_mistuning * randn(data.N, 1);
    
    k_aero = 0.01 * w0^2; 
    
    a_real_true = zeros(data.N, 1);
    a_imag_true = zeros(data.N, 1);
    
    a_real_true(end)     =  0.00 * k_aero; 
    a_real_true(2)     =  0.40 * k_aero; 
    a_real_true(1)   =  0.40 * k_aero; 
    a_real_true(3)     = -0.10 * k_aero; 
    a_real_true(end-1) = -0.10 * k_aero;
    
    a_imag_true(1)     =  1.00 * k_aero; 
    a_imag_true(2)     =  0.15 * k_aero; 
    a_imag_true(end)   =  0.15 * k_aero; 
    a_imag_true(3)     =  0.05 * k_aero; 
    a_imag_true(end-1) =  0.05 * k_aero; 

    % C. Forzante (fr, fi)
    % Onda viaggiante EO=6
    f_scalar_amp = 50.0; % Ampiezza fisica (arbitraria ma realistica per plot log)
    wave_shape = exp(1i * (0:data.N-1)' * (2*pi * engine_order / data.N));
    f_vec_complex = f_scalar_amp * wave_shape;
    
    f_real_true = real(f_vec_complex);
    f_imag_true = imag(f_vec_complex);
    
    % Assemblaggio theta_true [5N x 1]
    theta_true = [d0_true; a_real_true; a_imag_true; f_real_true; f_imag_true];
    
    %% 4. Calcolo Risposta (Forward Run)
    % Necessario calcolare la risposta "Mistuned" (Dati misurati) 
    % e la risposta "Tuned" (Riferimento ideale blu)
    
    % Helper per matrici circolanti (uguale a quella nel solver)
    mk_circ = @(v) toeplitz(v(end:-1:1), [v(end); v(1:end-1)]);
    
    A_aero_mat = mk_circ(a_real_true) + 1i * mk_circ(a_imag_true);
    D0_mat     = diag(d0_true);
    
    idx = data.idx_blades;
    
    % --- Sistema TUNED (Riferimento) ---
    A_tuned = data.K_nom;
    % Tuned ha d0=0, ma ha l'aerodinamica
    A_tuned(idx, idx) = A_tuned(idx, idx) + A_aero_mat; 
    
    % --- Sistema MISTUNED (Vero) ---
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
    
    %% 5. Aggiunta Rumore e Output
    
    % Parametro rumore per test_d0.m
    noise_percent = 0.00; % 5%
    max_amp = max(abs(y_mistuned_complex(:)));
    data.sigmaNoise = noise_percent * max_amp;
    
    % Generazione Rumore Complesso
    noise_mat = data.sigmaNoise * (randn(size(y_mistuned_complex)) + 1i*randn(size(y_mistuned_complex)));
    
    % data.y deve essere COMPLESSO per rotor_objective_with_grad
    % (Il residuo rt = data.y - xt deve avere senso vettoriale)
    data.y = y_mistuned_complex + noise_mat;
    
    % Salvataggio dati extra per il plotting in test_d0.m
    % test_d0 userà abs(data.y) per le linee grigie grazie al comando log10(abs(data.y))
    data.y_tuned = y_tuned_complex; % Risposta tuned (linea blu)
    data.y_clean = y_mistuned_complex; % Risposta vera senza rumore (linea nera)
    
    % Inizializzazione Dummy per data.scale (verrà sovrascritta da test_d0.m)
    % Serve per evitare errori se si passa data a funzioni prima di test_d0
    data.scale.d0 = 1; data.scale.ar = 1; data.scale.ai = 1; data.scale.f = 1;
    data.J_normalization = 1;
    
    fprintf('Dati generati. Sigma Noise: %g\n', data.sigmaNoise);
end
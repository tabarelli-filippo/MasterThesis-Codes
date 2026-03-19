function [J, grad] = rotor_objective_with_grad(theta_norm, data)
%ROTOR_OBJECTIVE_WITH_GRAD Evalutates J e Gradient 
    N = data.N;
    idx = data.idx_blades;
    
    % De-scaling
    s_d0 = data.scale.d0;
    s_ar = data.scale.ar;
    s_ai = data.scale.ai;
    s_f  = data.scale.f;
    
    d0_phys     = theta_norm(1:N) * s_d0;
    a_real_phys = theta_norm(N+1 : 2*N) * s_ar;
    a_imag_phys = theta_norm(2*N+1 : 3*N) * s_ai;
    f_real_phys = theta_norm(3*N+1 : 4*N) * s_f;
    f_imag_phys = theta_norm(4*N+1 : 5*N) * s_f;
    
    % Physical Matrices
    A_mat = data.K_nom; 
    M_mat = data.M_mat;
    N_sys = size(A_mat, 1);
    
    D0_blade = diag(d0_phys);
    A_aero   = make_circulant(a_real_phys) + 1i * make_circulant(a_imag_phys);
    A_mat(idx, idx) = A_mat(idx, idx) + D0_blade + A_aero;
    
    f_vec_full = zeros(N_sys, 1);
    f_vec_blade = f_real_phys + 1i * f_imag_phys;
    f_vec_full(idx) = f_vec_blade;
    
    % Solver
    [X, Lambda_mat] = eig(A_mat, M_mat);
    Lambda_diag = diag(Lambda_mat);
    Y = (M_mat * X) \ eye(N_sys); % Left Eigevectors normalized
    
    J_phys = 0; 
    GA_accum_full = zeros(N_sys, N_sys); 
    gf_accum_full = zeros(N_sys, 1); 
    
    % Residual
    for t = 1:length(data.omega)
        w_t = data.omega(t);
        inv_Lambda = 1 ./ (Lambda_diag - w_t^2);
        
        Yf = Y * f_vec_full;               
        xt_full = X * (inv_Lambda .* Yf);  
        xt_blades = xt_full(idx);
        
        rt = data.y(:, t) - xt_blades; 
        J_phys = J_phys + (rt' * rt); 

        % residual rt
        rt_full = zeros(N_sys, 1);
        rt_full(idx) = rt;
        
        % Fast computing
        step1 = rt_full' * X;         % r^H * X (H = 1)
        step2 = step1 .* (inv_Lambda.');  % r^H*X * (Lambda - omegaT^2*I)^(-1)
        row_vec_adj = step2 * Y;      % r^H*X*(Lambda - omegaT^2*I)^(-1) * Y
        
        % stock GA = sum( x_t * lambda_t^H )
        GA_accum_full = GA_accum_full + (xt_full * row_vec_adj);
        
        % stock gf = - sum( lambda_t^H )^T
        gf_accum_full = gf_accum_full - row_vec_adj.';
    end
    
    % Scaling OF, for convergence
    COST_SCALING = data.J_normalization;
    J = J_phys * COST_SCALING;

    % Gradient
    GA_blades = GA_accum_full(idx, idx);
    gf_blades = gf_accum_full(idx);

    grad_d0_phys = 2 * real(diag(GA_blades));
    grad_fr_phys =  2 * real(gf_blades);
    grad_fi_phys = -2 * imag(gf_blades);

    grad_ar_phys = zeros(N, 1);
    grad_ai_phys = zeros(N, 1);

    for j = 1:N
        p = 1:N;
        q = mod(p - j - 1, N) + 1;

        s = sum(GA_blades(sub2ind([N N], p, q)));
        grad_ai_phys(j) = -2 * imag(s);
    end

    for j = 1:N-1
        p = 1:N;
        q = mod(p - j - 1, N) + 1;

        s = sum(GA_blades(sub2ind([N N], p, q)));
        grad_ar_phys(j) = 2 * real(s);
    end
    
    grad_d0_norm = grad_d0_phys * s_d0;
    grad_ar_norm = grad_ar_phys * s_ar;
    grad_ai_norm = grad_ai_phys * s_ai;
    grad_fr_norm = grad_fr_phys * s_f;
    grad_fi_norm = grad_fi_phys * s_f;

    grad = [grad_d0_norm; grad_ar_norm; grad_ai_norm; grad_fr_norm; grad_fi_norm] * COST_SCALING;
end

%% Helper Functions
function C = make_circulant(v)
    C = toeplitz(v(end:-1:1), [v(end); v(1:end-1)]);
end

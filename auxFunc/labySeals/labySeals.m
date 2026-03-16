function [mdot_leakage,Kseal,Cseal] = labySeals(sealType,P_reserv,P_sump,Nt,Cr,L,T,Rs,B,U_inlet,omega,nu)
%LABYSEALS: computes leakage massflow, stiffness matrix and damping matrix
%   for Seals. Based on Childs-Scharrer (1896)
%
%INPUT:
% Nt: teeth number
% Cr: radial clearence
% L: pitch seal strip
% T: temperature
% B: tooth depth
% Cr: tooth

%% Fixed values
R = 287; %[J/Kg/K] gas
gamma = 1.4; % Cp/Cv gas
Nc = Nt-1; % Number of cavities

%% 1. PRESSURE DISTRIBUTION AND MASS LEAKAGE

M_trial = 1; % Mach number
pr_crit = pressureRatio(M_trial,gamma);
P_last_crit = P_sump * pr_crit;

alpha = 1 - (1 + 16.6*Cr/L)^(-2);
mu2 = (Nt/((1-alpha)*Nt + alpha))^(0.5);

mdot_limit = 0.510 * mu2 * P_last_crit * Cr / sqrt(R*T);

% P_reservoir needed for choking
[P_res_limit_calc, ~] = backwardPressureSolver(P_last_crit,Nt,mdot_limit,Cr,R,T,mu2,gamma);

% Solver
if P_reserv >= P_res_limit_calc
    % >> CHOKED FLOW
    obj_fun = @(P_guess) wrapperChokedError(P_guess, P_reserv,P_sump, Nt, Cr, R, T, mu2, gamma);

    options = optimset('Display','off');
    P_last_sol = fzero(obj_fun, P_last_crit * 1.1, options);

    [~, P_dist, mdot] = wrapperChokedError(P_last_sol, P_reserv,P_sump, Nt, Cr, R, T, mu2, gamma);

else
    % >> UNCHOKED FLOW
    obj_fun = @(P_guess) wrapperUnchokedError(P_guess, P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma);
    
    options = optimset('Display','off');

    term_drop = (P_reserv^2 - P_sump^2) / Nt;
    P_guess_parabolic = sqrt(P_reserv^2 - term_drop);
    
    lb = max(P_sump * 1.01, P_guess_parabolic * 0.90);
    
    ub = min(P_reserv * 0.99, P_guess_parabolic * 1.10);
    
    if lb >= ub
        lb = P_sump * 1.01;
        ub = P_reserv * 0.99;
    end
    
    try
        P_sol = fzero(obj_fun, [lb, ub], options);
    catch
        try
             lb_wide = P_sump * 1.05;
             ub_wide = P_reserv * 0.99;
             P_sol = fzero(obj_fun, [lb_wide, ub_wide], options);
        catch
             P_sol = fzero(obj_fun, P_guess_parabolic, options);
        end
    end
    [~, P_dist, mdot] = wrapperUnchokedError(P_sol, P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma);
end

mdot_leakage = mdot * 2*pi*Rs;

%% 2. TANGENTIAL VELOCITY
W_dist = zeros(1, Nt+1);
W_dist(1) = U_inlet;

% Hydraulic Diameter
Dh = 2 * (B + Cr) * L / (B + Cr + L);

if strcmpi(sealType, 'TOS')
    ar = 1; as = (L + 2*B) / L; % Stator Teeth
elseif strcmpi(sealType, 'TOR')
    ar = (L + 2*B) / L; as = 1; % Rotor Teeth
else
    error('sealType must be "TOS" or "TOR"');
end

% Blasius coefficents
ns = 0.079; ms = -0.25;
nr = 0.079; mr = -0.25;

U_rotor = omega * Rs;

for i = 1:Nc
    % Density and pressure cavity (i+1)
    P_cavity = P_dist(i+1);
    rho = P_cavity / (R * T);
    
    W_prev = W_dist(i);
    
    mom_eq = @(U) solveMomentum(U, W_prev, mdot, L, ar, as, rho, omega, ...
                                Rs, Cr, Dh, nu, ns, ms, nr, mr);
    
    lb = -0.1 * abs(U_rotor); 
    ub = 1.1 * abs(U_rotor);
    
    % if no pre swirl
    if abs(U_rotor) < 1e-3; lb = -10; ub = 10; end
    
    try
        % first try
        W_val = fzero(mom_eq, [lb, ub]);
    catch
        % Fallback
        try
            options = optimset('Display','on');
            W_val = fzero(mom_eq, W_prev,options); 
        catch
            W_val = 0.5 * U_rotor; 
        end
    end
    
    W_dist(i+1) = W_val;
end

W_dist(end) = W_dist(end-1);

%% 3. REDUCED MODEL COEFFICIENTS
[G_mat, X_mat] = sealCoeffs(P_dist, W_dist, mdot, Nc, Cr,...
    L, B, Rs, T, R, gamma, omega, nu, Dh, ar, as, ns, ms, nr, mr);

%% 4. GLOBAL MATRICES ASSEMBLING
Dim = 8 * Nc;
A_sys = zeros(Dim, Dim);
b_sys_a = zeros(Dim, 1);
b_sys_b = zeros(Dim, 1);
Ai0 = L * (Cr + B);

for j = 1:Nc
    % Ai-1 matrix
    P0i = P_dist(j+1); 
    V0i = W_dist(j+1);
    
    Aim1 = zeros(8,8);

    Aim1(1,2) = G_mat(j,4); Aim1(2,1) = G_mat(j,4); 
    Aim1(3,4) = G_mat(j,4); Aim1(4,3) = G_mat(j,4); 

    Aim1(5,2) = X_mat(j,4); Aim1(6,1) = X_mat(j,4);
    Aim1(7,4) = X_mat(j,4); Aim1(8,3) = X_mat(j,4);

    Aim1(5,6) = -mdot;      Aim1(6,5) = -mdot;
    Aim1(7,8) = -mdot;      Aim1(8,7) = -mdot;

    % Ai matrix
    Ai = zeros(8,8);
    Ai(1,1) = G_mat(j,1)*(omega + V0i/Rs);   Ai(2,2) = -G_mat(j,1)*(omega + V0i/Rs);
    Ai(3,3) = G_mat(j,1)*(-omega + V0i/Rs);  Ai(4,4) = -G_mat(j,1)*(-omega + V0i/Rs);

    Ai(1,2) = G_mat(j,3);   Ai(2,1) = G_mat(j,3);
    Ai(3,4) = G_mat(j,3);   Ai(4,3) = G_mat(j,3);

    Ai(5,2) = X_mat(j,3);   Ai(6,1) = X_mat(j,3);
    Ai(7,4) = X_mat(j,3);   Ai(8,3) = X_mat(j,3);

    Ai(5,1) = Ai0/Rs;   Ai(7,3) = Ai0/Rs; Ai(6,2) = -Ai0/Rs;    Ai(8,4) = -Ai0/Rs;

    Ai(5,5) = X_mat(j,1)*(omega + V0i/Rs);   Ai(6,6) = -X_mat(j,1)*(omega + V0i/Rs);
    Ai(7,7) = X_mat(j,1)*(-omega + V0i/Rs);  Ai(8,8) = -X_mat(j,1)*(-omega + V0i/Rs);

    Ai(5,6) = X_mat(j,2);   Ai(6,5) = X_mat(j,2);
    Ai(7,8) = X_mat(j,2);   Ai(8,7) = X_mat(j,2);

    Ai(1,5) = G_mat(j,1)*P0i/Rs;      Ai(3,7) = G_mat(j,1)*P0i/Rs;
    Ai(2,6) = -G_mat(j,1)*P0i/Rs;     Ai(4,8) = -G_mat(j,1)*P0i/Rs;

    % Ai-1 matrix
    Aip1 = zeros(8,8);
    Aip1(1,2) = G_mat(j,5); Aip1(2,1) = G_mat(j,5); 
    Aip1(3,4) = G_mat(j,5); Aip1(4,3) = G_mat(j,5);

    %B vector
    B_local = [G_mat(j,6)/2;  -G_mat(j,2)/2*(omega + V0i/Rs); G_mat(j,6)/2; G_mat(j,2)/2*(omega - V0i/Rs);
         -X_mat(j,5)/2;     0;      -X_mat(j,5)/2;     0; ];
    %C vector
    C_local = [-G_mat(j,6)/2;  G_mat(j,2)/2*(omega + V0i/Rs); G_mat(j,6)/2; G_mat(j,2)/2*(omega - V0i/Rs);
          X_mat(j,5)/2;     0;      -X_mat(j,5)/2;     0; ];
    
    % index
    r_idx = (j-1)*8 + (1:8);
    
    % Ai
    A_sys(r_idx, r_idx) = Ai;
    
    % if not first cavity
    if j > 1
        c_idx_prev = (j-2)*8 + (1:8);
        A_sys(r_idx, c_idx_prev) = Aim1;
    end
    
    % if not last cavity
    if j < Nc
        c_idx_next = (j)*8 + (1:8);
        A_sys(r_idx, c_idx_next) = Aip1;
    end

    b_sys_a(r_idx) = B_local;
    b_sys_b(r_idx) = C_local;
end

Z = A_sys \ [b_sys_a, b_sys_b];

% Initialization
F_a_sin = 0; F_a_cos = 0;
F_b_sin = 0; F_b_cos = 0;

for j = 1:Nc
    idx_base = (j-1)*8; 
    % Extract Pressure Perturbations from Column 1 (a perturbation)
    P_si_p_a = Z(idx_base + 1, 1);
    P_ci_p_a = Z(idx_base + 2, 1);
    P_si_m_a = Z(idx_base + 3, 1);
    P_ci_m_a = Z(idx_base + 4, 1);
    
    % Extract Pressure Perturbations from Column 2 (b perturbation)
    P_si_p_b = Z(idx_base + 1, 2);
    P_ci_p_b = Z(idx_base + 2, 2);
    P_si_m_b = Z(idx_base + 3, 2);
    P_ci_m_b = Z(idx_base + 4, 2);
    

    F_a_sin = F_a_sin + (P_si_p_a - P_si_m_a) * L;
    F_a_cos = F_a_cos + (P_ci_p_a + P_ci_m_a) * L;
    
    F_b_sin = F_b_sin + (P_si_p_b - P_si_m_b) * L;
    F_b_cos = F_b_cos + (P_ci_p_b + P_ci_m_b) * L;
end

K = pi * Rs * F_a_cos;
k = pi * Rs * F_b_sin;


C = - (pi * Rs / omega) * F_a_sin;
c = (pi * Rs / omega) * F_b_cos;

Kseal = [K, k; -k, K];
Cseal = [C, c; -c, C];
end


function r = pressureRatio(M,gamma)
r = (1 + (gamma-1)/2 * M^2) ^ (gamma/(gamma-1));
end



function [mdot_leakage, Kseal, Cseal] = labySeals(sealType, P_reserv, P_sump, Nt, Cr, L, T, Rs, B, U_inlet, omega, nu)
% LABYSEALS  Computes leakage mass flow rate and rotordynamic force
%   coefficients for a labyrinth seal using the Scharrer-Childs bulk-flow
%   model.
%
%   The solution proceeds in four stages:
%     1. Pressure distribution and leakage mass flow: the inter-tooth
%        cavity pressures are found by solving the continuity equation
%        across each tooth. The flow regime (choked or unchoked) is
%        determined automatically by comparing the reservoir pressure
%        against the choked-flow threshold. The Chaplygin discharge
%        coefficient is used for each tooth.
%     2. Tangential (swirl) velocity distribution: the circumferential
%        momentum equation is integrated cavity by cavity, accounting for
%        stator and rotor wall shear stresses via Blasius correlations.
%     3. Perturbation coefficient matrices: first-order perturbation of
%        the continuity and momentum equations yields the G (continuity)
%        and X (momentum) coefficient arrays for each cavity.
%     4. Global system assembly and solution: the coupled perturbation
%        equations for all cavities are assembled into a block-tridiagonal
%        linear system; the solution gives the pressure perturbations
%        induced by a harmonic rotor eccentricity. These are integrated
%        circumferentially to obtain the net rotordynamic forces, from
%        which K, k, C, c are extracted.
%
%   The gas is modelled as a perfect gas (R = 287 J/kg/K, gamma = 1.4).
%   Two seal geometries are supported:
%     'TOS' - Teeth On Stator
%     'TOR' - Teeth On Rotor
%
% SYNTAX
%   [mdot_leakage, Kseal, Cseal] = labySeals(sealType, P_reserv, P_sump, ...
%       Nt, Cr, L, T, Rs, B, U_inlet, omega, nu)
%
% INPUT ARGUMENTS
%   sealType  - (string) Seal geometry: 'TOS' (teeth on stator) or
%               'TOR' (teeth on rotor)
%   P_reserv  - (scalar double) Upstream reservoir pressure [Pa]
%   P_sump    - (scalar double) Downstream sump pressure [Pa]
%   Nt        - (scalar integer) Number of teeth (= number of pressure
%               drops). Number of cavities Nc = Nt - 1.
%   Cr        - (scalar double) Radial tooth clearance [m]
%   L         - (scalar double) Tooth pitch (axial spacing) [m]
%   T         - (scalar double) Gas temperature [K]
%   Rs        - (scalar double) Seal rotor radius [m]
%   B         - (scalar double) Tooth depth (cavity depth) [m]
%   U_inlet   - (scalar double) Inlet tangential (swirl) velocity [m/s].
%               Positive = in the direction of shaft rotation.
%   omega     - (scalar double) Shaft angular velocity [rad/s]
%   nu        - (scalar double) Gas kinematic viscosity [m²/s]
%
% OUTPUT ARGUMENTS
%   mdot_leakage - (scalar double) Leakage mass flow rate per unit
%                  circumference, integrated over 2*pi*Rs [kg/s]
%   Kseal        - (2x2 double) Rotordynamic stiffness matrix [N/m]:
%                    Kseal = [ K,  k]
%                            [-k,  K]
%                  K = direct stiffness; k = cross-coupled stiffness
%   Cseal        - (2x2 double) Rotordynamic damping matrix [N·s/m]:
%                    Cseal = [ C,  c]
%                            [-c,  C]
%                  C = direct damping; c = cross-coupled damping
%
% NOTES
%   - The Chaplygin discharge coefficient mu1 is evaluated per tooth from
%     the local pressure ratio using the Chaplygin approximation.
%   - Choked flow at a tooth is detected when the local pressure ratio
%     equals the critical ratio (M=1). The critical mdot is computed via
%     the isentropic choking formula: mdot_crit = 0.510*mu2*Pcrit*Cr/sqrt(R*T).
%   - Blasius coefficients for stator and rotor walls: ns = nr = 0.079,
%     ms = mr = -0.25 (turbulent pipe flow approximation).
%   - The hydraulic diameter is Dh = 2*(B+Cr)*L / (B+Cr+L).
%   - The perturbation system is linear in the eccentricity; K, k, C, c
%     are therefore independent of eccentricity amplitude.
%
% REFERENCES
%   Childs, D.W. & Scharrer, J.K. (1986). An Iwatsubo-based solution for
%     labyrinth seals: comparison to experimental results.
%     ASME Journal of Engineering for Gas Turbines and Power, 108(2), 325-331.
%
% EXAMPLE
%   [mdot, Ks, Cs] = labySeals('TOS', 7.6e5, 0.943e5, 16, ...
%       4.06e-4, 3.175e-3, 300, 0.0765, 3.175e-3, 0, 3000*pi/30, 0.144e-4);
%   fprintf('Leakage: %.4f kg/s\n', mdot);
%   fprintf('Direct stiffness K = %.2f N/m\n', Ks(1,1));
%
% SEE ALSO
%   bearingMatrix (type 12), backwardPressureSolver, forwardPressureSolver,
%   wrapperChokedError, wrapperUnchokedError, sealCoeffs, solveMomentum,
%   massFlowUnchoked, testSeals

%% Fixed thermodynamic properties (perfect gas model)
R     = 287;   % specific gas constant [J/kg/K]
gamma = 1.4;   % heat capacity ratio Cp/Cv
Nc    = Nt-1;  % number of inter-tooth cavities

%% Stage 1: Pressure distribution and leakage mass flow

M_trial  = 1;
pr_crit  = pressureRatio(M_trial, gamma);
P_last_crit = P_sump * pr_crit;

% contraction coefficient (Egli formula for tooth series)
alpha = 1 - (1 + 16.6*Cr/L)^(-2);
mu2   = (Nt / ((1-alpha)*Nt + alpha))^0.5;

% critical (choked) mass flow at the last tooth
mdot_limit = 0.510 * mu2 * P_last_crit * Cr / sqrt(R*T);

% upstream pressure needed to produce choked flow
[P_res_limit_calc, ~] = backwardPressureSolver(P_last_crit, Nt, mdot_limit, Cr, R, T, mu2, gamma);

% select choked or unchoked solver
if P_reserv >= P_res_limit_calc
    % --- CHOKED FLOW ---
    obj_fun = @(P_guess) wrapperChokedError(P_guess, P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma);
    options = optimset('Display','off');
    P_last_sol = fzero(obj_fun, P_last_crit * 1.1, options);
    [~, P_dist, mdot] = wrapperChokedError(P_last_sol, P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma);

else
    % --- UNCHOKED FLOW ---
    obj_fun = @(P_guess) wrapperUnchokedError(P_guess, P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma);
    options = optimset('Display','off');

    term_drop       = (P_reserv^2 - P_sump^2) / Nt;
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
            P_sol   = fzero(obj_fun, [lb_wide, ub_wide], options);
        catch
            P_sol   = fzero(obj_fun, P_guess_parabolic, options);
        end
    end
    [~, P_dist, mdot] = wrapperUnchokedError(P_sol, P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma);
end

mdot_leakage = mdot * 2*pi*Rs;

%% Stage 2: Tangential (swirl) velocity distribution

W_dist    = zeros(1, Nt+1);
W_dist(1) = U_inlet;

% hydraulic diameter of the inter-tooth cavity
Dh = 2 * (B + Cr) * L / (B + Cr + L);

% geometry factor for stator/rotor surface area ratios
if strcmpi(sealType, 'TOS')
    ar = 1;              % rotor wetted perimeter ratio
    as = (L + 2*B) / L;  % stator wetted perimeter ratio
elseif strcmpi(sealType, 'TOR')
    ar = (L + 2*B) / L;
    as = 1;
else
    error('sealType must be "TOS" or "TOR"');
end

% Blasius friction law coefficients (turbulent)
ns = 0.079; ms = -0.25;   % stator
nr = 0.079; mr = -0.25;   % rotor

U_rotor = omega * Rs;

for i = 1:Nc
    P_cavity = P_dist(i+1);
    rho      = P_cavity / (R * T);
    W_prev   = W_dist(i);

    mom_eq = @(U) solveMomentum(U, W_prev, mdot, L, ar, as, rho, omega, ...
                                Rs, Cr, Dh, nu, ns, ms, nr, mr);

    lb = -0.1 * abs(U_rotor);
    ub =  1.1 * abs(U_rotor);
    if abs(U_rotor) < 1e-3; lb = -10; ub = 10; end

    try
        W_val = fzero(mom_eq, [lb, ub]);
    catch
        try
            opt_fb = optimset('Display','on');
            W_val  = fzero(mom_eq, W_prev, opt_fb);
        catch
            W_val  = 0.5 * U_rotor;
        end
    end
    W_dist(i+1) = W_val;
end

W_dist(end) = W_dist(end-1);

%% Stage 3: Perturbation coefficients (G and X matrices)
[G_mat, X_mat] = sealCoeffs(P_dist, W_dist, mdot, Nc, Cr, ...
    L, B, Rs, T, R, gamma, omega, nu, Dh, ar, as, ns, ms, nr, mr);

%% Stage 4: Global system assembly and solution

Dim   = 8 * Nc;
A_sys = zeros(Dim, Dim);
b_sys_a = zeros(Dim, 1);   % forcing column: x-eccentricity (cos)
b_sys_b = zeros(Dim, 1);   % forcing column: y-eccentricity (sin)
Ai0 = L * (Cr + B);        % cavity transverse cross-section area

for j = 1:Nc
    P0i = P_dist(j+1);
    V0i = W_dist(j+1);

    % Upstream coupling block (i-1 → i)
    Aim1 = zeros(8,8);
    Aim1(1,2) = G_mat(j,4); Aim1(2,1) = G_mat(j,4);
    Aim1(3,4) = G_mat(j,4); Aim1(4,3) = G_mat(j,4);
    Aim1(5,2) = X_mat(j,4); Aim1(6,1) = X_mat(j,4);
    Aim1(7,4) = X_mat(j,4); Aim1(8,3) = X_mat(j,4);
    Aim1(5,6) = -mdot;      Aim1(6,5) = -mdot;
    Aim1(7,8) = -mdot;      Aim1(8,7) = -mdot;

    % Diagonal block (cavity i)
    Ai = zeros(8,8);
    Ai(1,1) =  G_mat(j,1)*(omega + V0i/Rs);
    Ai(2,2) = -G_mat(j,1)*(omega + V0i/Rs);
    Ai(3,3) =  G_mat(j,1)*(-omega + V0i/Rs);
    Ai(4,4) = -G_mat(j,1)*(-omega + V0i/Rs);
    Ai(1,2) = G_mat(j,3); Ai(2,1) = G_mat(j,3);
    Ai(3,4) = G_mat(j,3); Ai(4,3) = G_mat(j,3);
    Ai(5,2) = X_mat(j,3); Ai(6,1) = X_mat(j,3);
    Ai(7,4) = X_mat(j,3); Ai(8,3) = X_mat(j,3);
    Ai(5,1) = Ai0/Rs;  Ai(7,3) =  Ai0/Rs;
    Ai(6,2) = -Ai0/Rs; Ai(8,4) = -Ai0/Rs;
    Ai(5,5) =  X_mat(j,1)*(omega + V0i/Rs);
    Ai(6,6) = -X_mat(j,1)*(omega + V0i/Rs);
    Ai(7,7) =  X_mat(j,1)*(-omega + V0i/Rs);
    Ai(8,8) = -X_mat(j,1)*(-omega + V0i/Rs);
    Ai(5,6) = X_mat(j,2); Ai(6,5) = X_mat(j,2);
    Ai(7,8) = X_mat(j,2); Ai(8,7) = X_mat(j,2);
    Ai(1,5) =  G_mat(j,1)*P0i/Rs; Ai(3,7) =  G_mat(j,1)*P0i/Rs;
    Ai(2,6) = -G_mat(j,1)*P0i/Rs; Ai(4,8) = -G_mat(j,1)*P0i/Rs;

    % Downstream coupling block (i → i+1)
    Aip1 = zeros(8,8);
    Aip1(1,2) = G_mat(j,5); Aip1(2,1) = G_mat(j,5);
    Aip1(3,4) = G_mat(j,5); Aip1(4,3) = G_mat(j,5);

    % Local forcing vectors (a and b eccentricity perturbations)
    B_local = [G_mat(j,6)/2;  -G_mat(j,2)/2*(omega + V0i/Rs);
               G_mat(j,6)/2;   G_mat(j,2)/2*(omega - V0i/Rs);
              -X_mat(j,5)/2;   0; -X_mat(j,5)/2; 0];
    C_local = [-G_mat(j,6)/2;  G_mat(j,2)/2*(omega + V0i/Rs);
                G_mat(j,6)/2;  G_mat(j,2)/2*(omega - V0i/Rs);
                X_mat(j,5)/2;  0; -X_mat(j,5)/2; 0];

    r_idx = (j-1)*8 + (1:8);

    A_sys(r_idx, r_idx) = Ai;
    if j > 1
        c_idx_prev = (j-2)*8 + (1:8);
        A_sys(r_idx, c_idx_prev) = Aim1;
    end
    if j < Nc
        c_idx_next = j*8 + (1:8);
        A_sys(r_idx, c_idx_next) = Aip1;
    end

    b_sys_a(r_idx) = B_local;
    b_sys_b(r_idx) = C_local;
end

Z = A_sys \ [b_sys_a, b_sys_b];

% integrate pressure perturbations over the seal length → net force
F_a_sin = 0; F_a_cos = 0;
F_b_sin = 0; F_b_cos = 0;

for j = 1:Nc
    idx_base = (j-1)*8;
    P_si_p_a = Z(idx_base + 1, 1); P_ci_p_a = Z(idx_base + 2, 1);
    P_si_m_a = Z(idx_base + 3, 1); P_ci_m_a = Z(idx_base + 4, 1);
    P_si_p_b = Z(idx_base + 1, 2); P_ci_p_b = Z(idx_base + 2, 2);
    P_si_m_b = Z(idx_base + 3, 2); P_ci_m_b = Z(idx_base + 4, 2);

    F_a_sin = F_a_sin + (P_si_p_a - P_si_m_a) * L;
    F_a_cos = F_a_cos + (P_ci_p_a + P_ci_m_a) * L;
    F_b_sin = F_b_sin + (P_si_p_b - P_si_m_b) * L;
    F_b_cos = F_b_cos + (P_ci_p_b + P_ci_m_b) * L;
end

% extract rotordynamic coefficients
K =  pi * Rs * F_a_cos;
k =  pi * Rs * F_b_sin;
C = -(pi * Rs / omega) * F_a_sin;
c =  (pi * Rs / omega) * F_b_cos;

Kseal = [ K,  k; -k,  K];
Cseal = [ C,  c; -c,  C];
end

% -------------------------------------------------------------------------
function r = pressureRatio(M, gamma)
% PRESSURERATIO  Isentropic total-to-static pressure ratio at Mach number M.
r = (1 + (gamma-1)/2 * M^2) ^ (gamma/(gamma-1));
end

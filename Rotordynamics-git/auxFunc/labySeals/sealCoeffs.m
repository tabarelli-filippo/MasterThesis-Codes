function [G, X] = sealCoeffs(P_dist, U_dist, mdot, Nc, Cr, L, B, Rs, T, R, gamma, omega, nu, Dh, ar, as, ns, ms, nr, mr)
% SEALCOEFFS  Computes the first-order perturbation coefficient matrices
%   G (continuity) and X (momentum) for each inter-tooth cavity of a
%   labyrinth seal, as required by the Scharrer-Childs bulk-flow model.
%
%   For each cavity i (i = 1 to Nc), the linearised perturbation of the
%   steady-state continuity and circumferential momentum equations about
%   the zeroth-order (axisymmetric) solution yields six coefficients stored
%   in G(i,:) and X(i,:). These coefficients couple the perturbation
%   pressures and tangential velocities in adjacent cavities and are used
%   to assemble the global block-tridiagonal system in labySeals.
%
%   The Chaplygin variable discharge coefficient mu1 and its derivative
%   are evaluated at each tooth using the local zeroth-order pressure ratio.
%   Wall shear stresses on stator and rotor surfaces are computed via the
%   Blasius turbulent friction law.
%
% SYNTAX
%   [G, X] = sealCoeffs(P_dist, U_dist, mdot, Nc, Cr, L, B, Rs, T, R, ...
%       gamma, omega, nu, Dh, ar, as, ns, ms, nr, mr)
%
% INPUT ARGUMENTS
%   P_dist  - (1 x Nc+2 double) Zeroth-order pressure distribution [Pa].
%             P_dist(1) = P_reservoir, P_dist(end) = P_sump.
%             Cavity pressures are P_dist(2:Nc+1).
%   U_dist  - (1 x Nc+2 double) Zeroth-order tangential velocity [m/s].
%             U_dist(1) = inlet swirl; U_dist(2:Nc+1) = cavity values.
%   mdot    - (scalar double) Zeroth-order mass flow per unit circumference
%             [kg/(m·s)]
%   Nc      - (scalar integer) Number of inter-tooth cavities = Nt - 1
%   Cr      - (scalar double) Radial tooth clearance [m]
%   L       - (scalar double) Tooth pitch [m]
%   B       - (scalar double) Tooth depth (cavity depth) [m]
%   Rs      - (scalar double) Seal radius [m]
%   T       - (scalar double) Gas temperature [K]
%   R       - (scalar double) Gas specific constant [J/kg/K]
%   gamma   - (scalar double) Heat capacity ratio [-]
%   omega   - (scalar double) Shaft angular velocity [rad/s]
%   nu      - (scalar double) Gas kinematic viscosity [m²/s]
%   Dh      - (scalar double) Cavity hydraulic diameter [m]
%             = 2*(B+Cr)*L / (B+Cr+L)
%   ar      - (scalar double) Rotor wetted perimeter ratio [-]
%   as      - (scalar double) Stator wetted perimeter ratio [-]
%   ns, ms  - (scalar double) Blasius coefficients for stator friction
%             (ns = 0.079, ms = -0.25 for turbulent pipe flow)
%   nr, mr  - (scalar double) Blasius coefficients for rotor friction
%
% OUTPUT ARGUMENTS
%   G - (Nc x 6 double) Continuity equation perturbation coefficients.
%       Column meaning:
%         G(:,1): coefficient of pressure perturbation at cavity i (Ai0/R/T)
%         G(:,2): coefficient involving cavity pitch integral
%         G(:,3): coefficient of pressure perturbation from inlet tooth
%         G(:,4): coupling coefficient from upstream cavity (i-1)
%         G(:,5): coupling coefficient from downstream cavity (i+1)
%         G(:,6): direct forcing term (zero for this model)
%   X - (Nc x 6 double) Momentum equation perturbation coefficients.
%       Column meaning:
%         X(:,1): inertia term (Pi0*Ai0/R/T)
%         X(:,2): shear stress velocity derivative contribution
%         X(:,3): upstream-tooth pressure coupling
%         X(:,4): upstream-cavity pressure coupling
%         X(:,5): clearance geometry derivative term
%         X(:,6): (unused placeholder)
%
% REFERENCES
%   Childs, D.W. & Scharrer, J.K. (1986). An Iwatsubo-based solution for
%     labyrinth seals: comparison to experimental results.
%     ASME Journal of Engineering for Gas Turbines and Power, 108(2), 325-331.
%
% SEE ALSO
%   labySeals, solveMomentum

G = zeros(Nc, 6);
X = zeros(Nc, 6);

Ai0      = L * (Cr + B);       % cavity transverse cross-section [m²]
fact_gamma = (gamma-1)/gamma;  % isentropic exponent shorthand

for i = 1:Nc
    % zeroth-order pressures (cavity i and neighbours)
    Pi0   = P_dist(i+1);
    Pi_m1 = P_dist(i);
    Pi_p1 = P_dist(i+2);

    % zeroth-order tangential velocities
    Ui0   = U_dist(i+1);
    Ui_m1 = U_dist(i);

    %% Discharge coefficients and their derivatives

    % inlet tooth (i-1 → i)
    pr_in      = Pi_m1 / Pi0;
    beta_i     = pr_in^fact_gamma - 1;
    mu_i0      = pi / (pi + 2 - 5*beta_i + 2*beta_i^2);
    fact_beta_i = (-5 + 4*beta_i)/pi;

    % exit tooth (i → i+1)
    pr_out       = Pi0 / Pi_p1;
    beta_ip1     = pr_out^fact_gamma - 1;
    mu_ip10      = pi / (pi + 2 - 5*beta_ip1 + 2*beta_ip1^2);
    fact_beta_ip1 = (5 - 4*beta_ip1)/pi;

    %% Wall shear stresses (Blasius turbulent law)
    rho = Pi0 / (R * T);

    % stator shear stress
    re_s   = abs(Ui0) * Dh / nu;
    tau_s0 = 0.5 * rho * Ui0^2 * ns * re_s^ms * sign(Ui0);

    % rotor shear stress (velocity relative to rotor surface)
    V_rel  = omega * Rs - Ui0;
    re_r   = abs(V_rel) * Dh / nu;
    tau_r0 = 0.5 * rho * V_rel^2 * nr * re_r^mr * sign(V_rel);

    %% Continuity coefficients (G matrix)

    G(i,1) = Ai0 / (R * T);

    G(i,2) = Pi0 * L / (R * T);

    term_G3_1 =  Pi0 / (Pi0^2 - Pi_m1^2);
    term_G3_2 =  mu_ip10 * (-fact_beta_ip1) * fact_gamma * (1/Pi_p1) * pr_out^(1/gamma);
    term_G3_3 =  mu_i0   * fact_beta_i      * fact_gamma * (1/Pi0)  * pr_in^fact_gamma;
    G(i,3) = mdot * (term_G3_1 + term_G3_2 + term_G3_3);

    term_G4_1 = -Pi_m1 / (Pi_m1^2 - Pi0^2);
    term_G4_2 =  mu_i0 * fact_beta_i * (fact_gamma * (1/Pi0) * pr_in^(-1/gamma));
    G(i,4) = mdot * (term_G4_1 + term_G4_2);

    term_G5_1 = -Pi_p1 / (Pi0^2 - Pi_p1^2);
    term_G5_2 =  mu_ip10 * fact_beta_ip1 * (fact_gamma * (1/Pi_p1) * pr_out^fact_gamma);
    G(i,5) = -mdot * (term_G5_1 + term_G5_2);

    G(i,6) = 0;

    %% Momentum coefficients (X matrix)

    X(i,1) = Pi0 * Ai0 / (R * T);

    denom_s = Ui0;
    denom_r = omega*Rs - Ui0;
    X(i,2) = mdot + ...
        ((2 + ms)/denom_s) * L * as * tau_s0 + ...
        ((2 + mr)/denom_r) * L * ar * tau_r0;

    dU = Ui0 - Ui_m1;
    term_X3_visc_s =  L * as * tau_s0 / Pi0;
    term_X3_visc_r = -L * ar * tau_r0 / Pi0;
    term_X3_3 = -mdot * Pi_m1 / (Pi_m1^2 - Pi0^2) * dU;
    term_X3_4 =  mdot * fact_beta_i * (fact_gamma * (1/Pi0) * pr_in^fact_gamma) * dU;
    X(i,3) = term_X3_visc_s + term_X3_visc_r + term_X3_3 + term_X3_4;

    term_X4_1 =  mdot * Pi_m1 / (Pi_m1^2 - Pi0^2) * dU;
    term_X4_2 = -mdot * fact_beta_i * (fact_gamma * (1/Pi0) * pr_in^(-1/gamma)) * dU;
    X(i,4) = term_X4_1 + term_X4_2;

    term_X5_1 = -mdot / Cr * dU;
    term_X5_s = -(ms * Dh * L * as * tau_s0) / (2 * (B+Cr)^2);
    term_X5_r =  (mr * Dh * L * ar * tau_r0) / (2 * (B+Cr)^2);
    X(i,5) = term_X5_1 + term_X5_s + term_X5_r;
end
end

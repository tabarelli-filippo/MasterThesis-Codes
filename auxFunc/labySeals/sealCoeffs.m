function [G, X] = sealCoeffs(P_dist, U_dist, mdot, Nc, Cr, L, B, Rs, T, R, gamma, omega, nu, Dh, ar, as, ns, ms, nr, mr)
%SEALCOEFFS: Perturbation Coefficients


G = zeros(Nc, 6);
X = zeros(Nc, 6);

% Cavity transverse area
Ai0 = L * (Cr + B);

lambda_visc = 0;

fact_gamma = (gamma-1)/gamma;

%% Loop over cavities (i = 1 to Nc)
for i = 1:Nc
    
    Pi0   = P_dist(i+1); % i th cavity pressure
    Pi_m1 = P_dist(i);   % i-1 th cavity pressure
    Pi_p1 = P_dist(i+2); % i+1 th cavity pressure

    Ui0   = U_dist(i+1); % i th tangential velocity 
    Ui_m1 = U_dist(i);   % i-1 th tangential velocity 

    % Auxiliary parameters
    pr_in = Pi_m1 / Pi0;
    beta_i = (pr_in)^fact_gamma - 1; 
    mu_i0 = pi / (pi + 2 - 5*beta_i + 2*beta_i^2); 

    % exit area
    pr_out = Pi0 / Pi_p1;
    beta_ip1 = (pr_out)^fact_gamma - 1;
    mu_ip10 = pi / (pi + 2 - 5*beta_ip1 + 2*beta_ip1^2);

    fact_beta_i   = (-5 + 4*beta_i)/pi;
    fact_beta_ip1 = (5 - 4*beta_ip1)/pi;

    %%  Shear Stress 
    rho = Pi0 / (R * T);
    % Stator
    re_s = (abs(Ui0) * Dh / nu);
    tau_s0 = 0.5 * rho * Ui0^2 * ns * (re_s)^ms * sign(Ui0);

    % Rotor
    V_rel = omega * Rs - Ui0;
    re_r = (abs(V_rel) * Dh / nu);
    tau_r0 = 0.5 * rho * V_rel^2 * nr * (re_r)^mr * sign(V_rel);

    %% Continuity coefficients

    G(i, 1) = Ai0 / (R * T);

    G(i, 2) = Pi0 * L / (R * T);

    term_G3_1 = Pi0 / (Pi0^2 - Pi_m1^2);
    term_G3_2 = mu_ip10 * (-fact_beta_ip1) * fact_gamma * (1/Pi_p1) * (pr_out)^(1/gamma);
    term_G3_3 = mu_i0 * fact_beta_i * fact_gamma * (1/Pi0) * (pr_in)^fact_gamma;
    G(i, 3) = mdot * (term_G3_1 + term_G3_2 + term_G3_3 );

    term_G4_1 = -Pi_m1 / (Pi_m1^2 - Pi0^2);
    term_G4_2 = mu_i0 * fact_beta_i * (fact_gamma * (1/Pi0) * (pr_in)^(-1/gamma));

    G(i, 4) = mdot * (term_G4_1 + term_G4_2);

    term_G5_1 = -Pi_p1 / (Pi0^2 - Pi_p1^2);
    term_G5_2 = mu_ip10 * fact_beta_ip1 * (fact_gamma * (1/Pi_p1) * (pr_out)^(fact_gamma));

    G(i, 5) = -mdot * (term_G5_1 + term_G5_2);
    G(i, 6) = 0;
    %% 4. Momentum coefficents 

    X(i, 1) = Pi0 * Ai0 / (R * T);

    denom_s = Ui0;
    denom_r = omega*Rs - Ui0;
    X(i, 2) = mdot + ...
        ((2 + ms)/denom_s) * L * as * tau_s0 + ...
        ((2 + mr)/denom_r) * L * ar * tau_r0;


    dU = (Ui0 - Ui_m1);
    term_X3_visc_s = L * as * tau_s0 / Pi0;
    term_X3_visc_r = - L * ar * tau_r0 / Pi0;
    term_X3_3 = - mdot * Pi_m1 / (Pi_m1^2 - Pi0^2) * dU;
    term_X3_4 =  mdot * fact_beta_i * (fact_gamma * (1/Pi0) * (pr_in)^(fact_gamma))*dU;

    X(i, 3) = term_X3_visc_s + term_X3_visc_r + term_X3_3 + term_X3_4;


    term_X4_1 = mdot * Pi_m1 / (Pi_m1^2 - Pi0^2) * dU;
    term_X4_2 = - mdot * fact_beta_i * (fact_gamma * (1/Pi0) * (pr_in)^(-1/gamma))*dU;

    X(i, 4) = term_X4_1 + term_X4_2;

    term_X5_1 = - mdot / Cr * dU;
    term_X5_s = - (ms * Dh * L * as * tau_s0) / (2 * (B+Cr)^2) ;
    term_X5_r = (mr * Dh * L * ar * tau_r0) / (2 * (B+Cr)^2) ;

    X(i, 5) = term_X5_1 + term_X5_s + term_X5_r;

end
end
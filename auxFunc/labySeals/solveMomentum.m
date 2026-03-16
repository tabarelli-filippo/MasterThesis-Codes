function F = solveMomentum(W, W_prev, mdot, L, ar, as, rho, omega, Rs,...
    Cr, Dh, nu, ns, ms, nr, mr)
%SOLVEMOMENTUM: auxiliary function for labySeals.m
% defines equation to be solved for W tangential velocity
    term_s = (abs(W) * Dh / nu)^ms;
    tau_s = 0.5 * rho * W^2 * ns * term_s * sign(W);
    
    % rotor shear stress
    V_rel = omega*Rs - W;
    term_r = (abs(V_rel) * Dh / nu)^mr;
    tau_r = 0.5*rho*V_rel^2 * nr * term_r * sign(V_rel);
    
    % residual
    F = mdot*(W-W_prev) - L*(ar*tau_r - as*tau_s);
end

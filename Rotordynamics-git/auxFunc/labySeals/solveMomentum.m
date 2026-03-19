function F = solveMomentum(W, W_prev, mdot, L, ar, as, rho, omega, Rs, Cr, Dh, nu, ns, ms, nr, mr)
% SOLVEMOMENTUM  Evaluates the residual of the steady-state circumferential
%   momentum equation for one inter-tooth cavity of a labyrinth seal.
%
%   For each cavity i, the tangential momentum balance between axial planes
%   at positions z_i-1 and z_i reads:
%
%       mdot * (W_i - W_{i-1}) = L * (ar*tau_r - as*tau_s)
%
%   where tau_r and tau_s are the rotor and stator wall shear stresses
%   computed via the Blasius turbulent friction law:
%       tau = 0.5 * rho * V² * n * (|V|*Dh/nu)^m * sign(V)
%
%   This function returns the residual F (= left-hand side minus
%   right-hand side), which is passed to fzero in labySeals to find the
%   cavity tangential velocity W_i that satisfies the momentum balance.
%
% SYNTAX
%   F = solveMomentum(W, W_prev, mdot, L, ar, as, rho, omega, Rs, ...
%       Cr, Dh, nu, ns, ms, nr, mr)
%
% INPUT ARGUMENTS
%   W      - (scalar double) Trial tangential velocity in cavity i [m/s].
%             This is the unknown being solved for by fzero.
%   W_prev - (scalar double) Tangential velocity in cavity i-1 [m/s]
%             (known from previous iteration or inlet condition)
%   mdot   - (scalar double) Mass flow per unit circumference [kg/(m·s)]
%   L      - (scalar double) Tooth pitch [m]
%   ar     - (scalar double) Rotor wetted perimeter ratio [-]
%   as     - (scalar double) Stator wetted perimeter ratio [-]
%   rho    - (scalar double) Gas density in cavity i [kg/m³]
%             = P_cavity / (R*T)
%   omega  - (scalar double) Shaft angular velocity [rad/s]
%   Rs     - (scalar double) Seal rotor radius [m]
%   Cr     - (scalar double) Radial clearance [m] (not used directly
%             here but available for extended models)
%   Dh     - (scalar double) Cavity hydraulic diameter [m]
%   nu     - (scalar double) Gas kinematic viscosity [m²/s]
%   ns, ms - (scalar double) Blasius stator friction coefficients
%            (ns = 0.079, ms = -0.25)
%   nr, mr - (scalar double) Blasius rotor friction coefficients
%
% OUTPUT ARGUMENTS
%   F - (scalar double) Momentum residual [N/m²]. Zero when W is the
%       correct tangential velocity satisfying the steady-state momentum
%       balance for cavity i.
%
% SEE ALSO
%   labySeals, sealCoeffs

    % stator shear stress (velocity W relative to stator at rest)
    term_s = (abs(W) * Dh / nu)^ms;
    tau_s  = 0.5 * rho * W^2 * ns * term_s * sign(W);

    % rotor shear stress (velocity relative to rotating rotor surface)
    V_rel  = omega*Rs - W;
    term_r = (abs(V_rel) * Dh / nu)^mr;
    tau_r  = 0.5 * rho * V_rel^2 * nr * term_r * sign(V_rel);

    % residual: momentum flux = net wall shear force
    F = mdot*(W - W_prev) - L*(ar*tau_r - as*tau_s);
end

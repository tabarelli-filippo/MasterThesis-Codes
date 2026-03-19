function [P_res_calc, P_distribution] = backwardPressureSolver(P_last_cavity, Nt, mdot, Cr, R, T, mu2, gamma)
% BACKWARDPRESSURESOLVER  Reconstructs the inter-tooth pressure distribution
%   of a labyrinth seal by marching upstream from a known cavity pressure.
%
%   Given the pressure in the last inter-tooth cavity (immediately upstream
%   of the sump), the function iterates backward through each tooth to
%   compute the pressure in each preceding cavity. The Chaplygin discharge
%   coefficient is applied at each tooth, and fzero is used to invert the
%   mass-flow equation for the upstream pressure.
%
%   This function is an auxiliary solver called by wrapperChokedError and
%   internally within labySeals.
%
% SYNTAX
%   [P_res_calc, P_distribution] = backwardPressureSolver(P_last_cavity, ...
%       Nt, mdot, Cr, R, T, mu2, gamma)
%
% INPUT ARGUMENTS
%   P_last_cavity - (scalar double) Pressure in the last inter-tooth
%                   cavity (immediately before the sump tooth) [Pa]
%   Nt            - (scalar integer) Number of teeth traversed in this
%                   backward pass (= total teeth - 1 for the choked case)
%   mdot          - (scalar double) Mass flow rate per unit circumference
%                   [kg/(m·s)]
%   Cr            - (scalar double) Radial tooth clearance [m]
%   R             - (scalar double) Gas specific constant [J/kg/K]
%   T             - (scalar double) Gas temperature [K]
%   mu2           - (scalar double) Contraction coefficient (Egli series
%                   formula), evaluated in labySeals [-]
%   gamma         - (scalar double) Heat capacity ratio Cp/Cv [-]
%
% OUTPUT ARGUMENTS
%   P_res_calc    - (scalar double) Computed reservoir pressure [Pa]
%                   after marching back through all Nt teeth from
%                   P_last_cavity. Used as residual in wrapperChokedError.
%   P_distribution- (1 x Nt+1 double) Pressure vector [Pa] with
%                   P_distribution(1) = P_res_calc (reservoir side) and
%                   P_distribution(end) = P_last_cavity (sump side).
%
% ALGORITHM
%   For each cavity i (from Nt-1 down to 1), the upstream pressure Pin is
%   found by solving the implicit equation:
%       (mdot / (mu1*mu2*Cr))^2 * R*T = Pin^2 - Pout^2
%   where mu1 = pi / (pi + 2 - 5*s + 2*s^2) is the Chaplygin discharge
%   coefficient with s = (Pin/Pout)^((gamma-1)/gamma) - 1.
%   The search bracket is [Pout*1.01, Pout*100].
%
% SEE ALSO
%   labySeals, wrapperChokedError, forwardPressureSolver

    P_dist_internal = zeros(1, Nt);
    P_dist_internal(Nt) = P_last_cavity;

    current_P_downstream = P_last_cavity;

    % march backward through cavities Nt-1 → 1
    for i = (Nt-1):-1:1
        P_upstream = pressure_backward(current_P_downstream, Cr, R, T, mu2, gamma, mdot);
        P_dist_internal(i) = P_upstream;
        current_P_downstream = P_upstream;
    end

    % one additional step to recover reservoir pressure
    P_res_calc   = pressure_backward(current_P_downstream, Cr, R, T, mu2, gamma, mdot);
    P_distribution = [P_res_calc, P_dist_internal];
end

% -------------------------------------------------------------------------
function Pin = pressure_backward(P_downstream, Cr, R, T, mu2, gamma, mdot)
% Invert the mass-flow equation to find upstream pressure given downstream.
    Pin_res = @(P_upstream) Pin_wrapperError(P_upstream, P_downstream, Cr, R, T, mu2, gamma, mdot);
    Pin = fzero(Pin_res, [P_downstream*1.01, P_downstream*100]);
end

% -------------------------------------------------------------------------
function res = Pin_wrapperError(Pin, Pout, Cr, R, T, mu2, gamma, mdot)
% Residual of the tooth mass-flow equation (backward direction).
% Returns a large value if Pin <= Pout to enforce a physical bracket.
    if Pin <= Pout; res = 1e20; return; end

    s    = (Pin/Pout)^((gamma-1)/gamma) - 1;
    mu1  = pi / (pi + 2 - 5*s + 2*s^2);   % Chaplygin coefficient

    kin_term       = (real(mdot)/mu1/mu2/Cr)^2 * R * T;
    delta_Pressure = Pin^2 - Pout^2;

    res = kin_term - delta_Pressure;
end

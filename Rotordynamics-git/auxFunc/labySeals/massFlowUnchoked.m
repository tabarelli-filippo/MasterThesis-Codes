function mdot = massFlowUnchoked(P_up, P_down, Cr, R, T, gamma, mu2)
% MASSFLOWUNCHOKED  Computes the unchoked leakage mass flow rate per unit
%   circumference through a single labyrinth seal tooth using the Chaplygin
%   compressible discharge model.
%
%   The mass flow is evaluated from the isentropic pressure-drop equation
%   with the Chaplygin variable-discharge coefficient mu1:
%
%       mdot = Cr * mu1 * mu2 * sqrt((P_up² - P_down²) / (R*T))
%
%   where mu1 is a function of the local pressure ratio:
%       s    = (P_up/P_down)^((gamma-1)/gamma) - 1
%       mu1  = pi / (pi + 2 - 5*s + 2*s^2)
%
%   This function is called by wrapperUnchokedError within labySeals to
%   compute the mass flow once the first-cavity pressure is known.
%
% SYNTAX
%   mdot = massFlowUnchoked(P_up, P_down, Cr, R, T, gamma, mu2)
%
% INPUT ARGUMENTS
%   P_up   - (scalar double) Upstream (reservoir-side) pressure [Pa]
%   P_down - (scalar double) Downstream (first-cavity) pressure [Pa]
%   Cr     - (scalar double) Radial tooth clearance [m]
%   R      - (scalar double) Gas specific constant [J/kg/K]
%   T      - (scalar double) Gas temperature [K]
%   gamma  - (scalar double) Heat capacity ratio Cp/Cv [-]
%   mu2    - (scalar double) Series contraction coefficient (Egli formula)
%            [-], evaluated once in labySeals for the whole seal
%
% OUTPUT ARGUMENTS
%   mdot - (scalar double) Mass flow rate per unit circumference
%          [kg/(m·s)]. Returns 0 if P_up <= P_down (no flow or
%          backflow condition).
%
% NOTES
%   - The formula is valid for unchoked (subsonic) flow; for choked flow
%     use the critical mass flow formula in labySeals directly.
%   - mu1 approaches 0.611 (Kirchhoff value) for small pressure drops
%     and decreases at higher pressure ratios.
%
% SEE ALSO
%   labySeals, wrapperUnchokedError, backwardPressureSolver

    if P_up <= P_down
        mdot = 0;
        return;
    end

    pr_val = P_down / P_up;
    s    = (1/pr_val)^((gamma-1)/gamma) - 1;
    mu1  = pi / (pi + 2 - 5*s + 2*s^2);    % Chaplygin coefficient

    mdot = Cr * mu1 * mu2 * sqrt((P_up^2 - P_down^2) / R / T);
end

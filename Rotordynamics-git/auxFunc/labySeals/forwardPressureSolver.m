function [P_end_calc, P_distribution] = forwardPressureSolver(P_start, Nt, mdot, Cr, R, T, mu2, gamma, P_targets)
% FORWARDPRESSURESOLVER  Reconstructs the inter-tooth pressure distribution
%   of a labyrinth seal by marching downstream from a known first-cavity
%   pressure.
%
%   Given the pressure immediately downstream of the reservoir tooth
%   (P_start), the function steps forward cavity by cavity to compute the
%   pressure at each subsequent inter-tooth location. At each step, fzero
%   inverts the Chaplygin-based mass-flow equation for the downstream
%   pressure, using P_targets as warm-start guesses for faster convergence.
%
%   This function is an auxiliary solver called by wrapperUnchokedError
%   within labySeals.
%
% SYNTAX
%   [P_end_calc, P_distribution] = forwardPressureSolver(P_start, Nt, ...
%       mdot, Cr, R, T, mu2, gamma, P_targets)
%
% INPUT ARGUMENTS
%   P_start    - (scalar double) Pressure in the first inter-tooth cavity,
%                immediately downstream of the reservoir tooth [Pa]
%   Nt         - (scalar integer) Number of teeth (forward steps) to
%                traverse
%   mdot       - (scalar double) Mass flow rate per unit circumference
%                [kg/(m·s)]
%   Cr         - (scalar double) Radial tooth clearance [m]
%   R          - (scalar double) Gas specific constant [J/kg/K]
%   T          - (scalar double) Gas temperature [K]
%   mu2        - (scalar double) Series contraction coefficient [-]
%   gamma      - (scalar double) Heat capacity ratio Cp/Cv [-]
%   P_targets  - (1 x Nt+1 double) Initial pressure guesses for each
%                cavity, used to warm-start fzero at each step [Pa].
%                Typically a parabolic (equal-drop in P²) distribution.
%
% OUTPUT ARGUMENTS
%   P_end_calc    - (scalar double) Computed pressure at the last cavity
%                   after traversing all Nt teeth [Pa]. Used as the
%                   residual in wrapperUnchokedError (should equal P_sump).
%   P_distribution- (1 x Nt+1 double) Full cavity pressure vector [Pa]:
%                   P_distribution(1) = P_start,
%                   P_distribution(end) = P_end_calc.
%
% ALGORITHM
%   At each tooth i (1 to Nt), the downstream pressure Pout is found by
%   solving:
%       (mdot / (mu1*mu2*Cr))^2 * R*T = Pin^2 - Pout^2
%   where mu1 is the Chaplygin coefficient evaluated at the local pressure
%   ratio. The search uses a warm-start from P_targets(i+1) when it lies
%   within the valid bracket [Pin*1e-4, Pin*0.99999].
%   A bracketed fzero search over the full range is used as fallback.
%
% SEE ALSO
%   labySeals, wrapperUnchokedError, backwardPressureSolver

    P_distribution    = zeros(1, Nt + 1);
    P_distribution(1) = P_start;

    P_upstream = P_start;

    for i = 1:Nt
        P_guess_local = P_targets(i+1);
        P_downstream  = pressure_forward(P_upstream, Cr, R, T, mu2, gamma, mdot, P_guess_local);
        P_distribution(i+1) = P_downstream;
        P_upstream = P_downstream;
    end

    P_end_calc = P_distribution(end);
end

% -------------------------------------------------------------------------
function Pout = pressure_forward(P_upstream, Cr, R, T, mu2, gamma, mdot, P_guess_local)
% Invert the mass-flow equation to find downstream pressure given upstream.
    Pout_res = @(P_down) Pin_wrapperError_Fwd(P_upstream, P_down, Cr, R, T, mu2, gamma, mdot);

    ub = P_upstream * 0.99999;
    lb = P_upstream * 1e-4;

    if P_guess_local < ub && P_guess_local > lb
        try
            Pout = fzero(Pout_res, P_guess_local);
            return;
        catch
            lb_search = max(lb, P_guess_local * 0.5);
            ub_search = min(ub, P_guess_local * 1.5);
            if lb_search >= ub_search; lb_search = lb; ub_search = ub; end
        end
    else
        lb_search = lb;
        ub_search = ub;
    end

    try
        Pout = fzero(Pout_res, [lb_search, ub_search]);
    catch
        Pout = lb;  % fallback: minimum physically admissible value
    end
end

% -------------------------------------------------------------------------
function res = Pin_wrapperError_Fwd(Pin, Pout, Cr, R, T, mu2, gamma, mdot)
% Residual of the tooth mass-flow equation (forward direction).
    if Pout >= Pin;  res =  1e20; return; end
    if Pout <= 0;    res = -1e20; return; end

    s    = (Pin/Pout)^((gamma-1)/gamma) - 1;
    mu1  = pi / (pi + 2 - 5*s + 2*s^2);

    kin_term       = (real(mdot)/mu1/mu2/Cr)^2 * R * T;
    delta_Pressure = Pin^2 - Pout^2;

    res = kin_term - delta_Pressure;
end

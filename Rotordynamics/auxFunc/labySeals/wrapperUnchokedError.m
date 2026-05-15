function [err, P_distribution, mdot] = wrapperUnchokedError(P_first_guess, P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma)
% WRAPPERUNCHOKEDERROR  Residual function for the unchoked-flow pressure
%   distribution solver in labySeals.
%
%   When the labyrinth seal operates in the unchoked regime, the mass flow
%   is determined by the pressure drop across the first tooth (from the
%   reservoir to the first inter-tooth cavity at P_first_guess) using the
%   Chaplygin unchoked mass-flow formula. Given this mdot, the pressure
%   distribution is propagated forward through all Nt-1 downstream teeth
%   via forwardPressureSolver. The residual is the difference between the
%   computed last-cavity pressure and the actual sump pressure.
%
%   This function is passed to fzero in labySeals to find the P_first_guess
%   that makes the forward-computed pressure match P_sump.
%
% SYNTAX
%   [err, P_distribution, mdot] = wrapperUnchokedError(P_first_guess, ...
%       P_reserv, P_sump, Nt, Cr, R, T, mu2, gamma)
%
% INPUT ARGUMENTS
%   P_first_guess - (scalar double) Trial pressure in the first inter-tooth
%                   cavity (immediately downstream of the reservoir tooth)
%                   [Pa]
%   P_reserv      - (scalar double) Reservoir (upstream) pressure [Pa]
%   P_sump        - (scalar double) Sump (downstream) pressure [Pa]
%   Nt            - (scalar integer) Total number of teeth
%   Cr            - (scalar double) Radial tooth clearance [m]
%   R             - (scalar double) Gas specific constant [J/kg/K]
%   T             - (scalar double) Gas temperature [K]
%   mu2           - (scalar double) Series contraction coefficient [-]
%   gamma         - (scalar double) Heat capacity ratio Cp/Cv [-]
%
% OUTPUT ARGUMENTS
%   err          - (scalar double) Residual: P_sump_computed - P_sump [Pa].
%                  Zero when P_first_guess produces the correct last-cavity
%                  pressure.
%   P_distribution - (1 x Nt+1 double) Complete pressure distribution [Pa]:
%                   P_distribution(1) = P_reserv,
%                   P_distribution(end) = P_sump_computed.
%   mdot         - (scalar double) Unchoked mass flow per unit circumference
%                  [kg/(m·s)] corresponding to P_first_guess.
%
% NOTES
%   - A parabolic (equal-drop in P²) distribution is used as the warm-start
%     guess vector for forwardPressureSolver.
%   - The reservoir pressure P_reserv is prepended to the partial
%     distribution returned by forwardPressureSolver.
%
% SEE ALSO
%   labySeals, forwardPressureSolver, wrapperChokedError, massFlowUnchoked

    % unchoked mass flow across the first tooth
    mdot = massFlowUnchoked(P_reserv, P_first_guess, Cr, R, T, gamma, mu2);

    % parabolic warm-start guesses for the forward solver (equal P² drops)
    P_targets    = zeros(1, Nt);
    P_targets(1) = P_first_guess;
    term_drop    = (P_first_guess^2 - P_sump^2) / (Nt-1);
    for k = 1:(Nt-1)
        P_targets(k+1) = sqrt(max(0, P_first_guess^2 - k * term_drop));
    end

    % forward pressure reconstruction through Nt-1 downstream teeth
    [P_sump_calc, P_dist_partial] = forwardPressureSolver(P_first_guess, Nt-1, mdot, Cr, R, T, mu2, gamma, P_targets);

    % prepend reservoir pressure to complete the distribution
    P_distribution = [P_reserv, P_dist_partial];

    err = P_sump_calc - P_sump;
end

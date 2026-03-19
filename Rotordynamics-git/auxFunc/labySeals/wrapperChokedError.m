function [err, P_dist, mdot] = wrapperChokedError(P_last_guess, P_reservoir, P_sump, Nt, Cr, R, T, mu2, gamma)
% WRAPPERCHOKEDERROR  Residual function for the choked-flow pressure
%   distribution solver in labySeals.
%
%   When the last tooth of the labyrinth seal is choked, the mass flow is
%   uniquely determined by the pressure immediately upstream of the last
%   tooth (P_last_guess) via the isentropic choking formula:
%       mdot = 0.510 * mu2 * P_last_guess * Cr / sqrt(R*T)
%
%   Given this mdot, the pressure distribution is reconstructed backward
%   from P_last_guess through Nt-1 upstream teeth using
%   backwardPressureSolver. The residual is the difference between the
%   computed reservoir pressure and the actual reservoir pressure P_reservoir.
%
%   This function is passed to fzero in labySeals to find the P_last_guess
%   that makes the residual zero, i.e., makes the computed upstream pressure
%   match P_reservoir.
%
% SYNTAX
%   [err, P_dist, mdot] = wrapperChokedError(P_last_guess, P_reservoir, ...
%       P_sump, Nt, Cr, R, T, mu2, gamma)
%
% INPUT ARGUMENTS
%   P_last_guess - (scalar double) Trial pressure immediately upstream of
%                  the last (choked) tooth [Pa]
%   P_reservoir  - (scalar double) Actual reservoir pressure [Pa]
%   P_sump       - (scalar double) Sump (downstream) pressure [Pa]
%   Nt           - (scalar integer) Total number of teeth
%   Cr           - (scalar double) Radial tooth clearance [m]
%   R            - (scalar double) Gas specific constant [J/kg/K]
%   T            - (scalar double) Gas temperature [K]
%   mu2          - (scalar double) Series contraction coefficient [-]
%   gamma        - (scalar double) Heat capacity ratio Cp/Cv [-]
%
% OUTPUT ARGUMENTS
%   err    - (scalar double) Residual: P_res_computed - P_reservoir [Pa].
%            Zero when P_last_guess yields the correct upstream pressure.
%   P_dist - (1 x Nt+1 double) Complete pressure distribution [Pa]:
%            P_dist(1) = P_reservoir (computed), P_dist(end) = P_sump.
%   mdot   - (scalar double) Choked mass flow rate per unit circumference
%            [kg/(m·s)] corresponding to P_last_guess.
%
% SEE ALSO
%   labySeals, backwardPressureSolver, wrapperUnchokedError

    % choked mass flow from last-tooth pressure (isentropic critical condition)
    mdot = 0.510 * mu2 * P_last_guess * Cr / sqrt(R*T);

    % backward pressure reconstruction through Nt-1 upstream teeth
    [P_res_calc, P_dist_partial] = backwardPressureSolver(P_last_guess, Nt-1, mdot, Cr, R, T, mu2, gamma);

    % append sump pressure to complete the distribution
    P_dist = [P_dist_partial, P_sump];

    err = P_res_calc - P_reservoir;
end

function [err, P_distribution, mdot] = wrapperUnchokedError(P_first_guess, P_reserv,P_sump, Nt, Cr, R, T, mu2, gamma)
% WRAPPERUNCHOCKEDERRROR: auxiliary function for labySeals.m
% defines equation to be solved for pressure distribution - unchoked
% massflow
    mdot = massFlowUnchoked(P_reserv, P_first_guess, Cr, R, T, gamma, mu2);

    idx = 0:Nt; 
    P_targets = zeros(1, Nt);
    P_targets(1) = P_first_guess;
    
    
    term_drop = (P_first_guess^2 - P_sump^2) / (Nt-1);
    
    for k = 1:(Nt-1)
        P_targets(k+1) = sqrt( max(0, P_first_guess^2 - k * term_drop) );
    end

    [P_sump_calc, P_dist_partial] = forwardPressureSolver(P_first_guess, Nt-1, mdot, Cr, R, T, mu2, gamma, P_targets);
    
    P_distribution = [P_reserv, P_dist_partial];
    err = P_sump_calc - P_sump;
end
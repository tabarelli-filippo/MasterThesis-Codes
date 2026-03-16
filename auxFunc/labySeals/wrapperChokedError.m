function [err, P_dist, mdot] = wrapperChokedError(P_last_guess, P_reservoir,P_sump, Nt, Cr, R, T, mu2, gamma)
% WRAPPERCHOCKEDERRROR: auxiliary function for labySeals.m
% defines equation to be solved for pressure distribution - choked massflow
    mdot = 0.510 * mu2 * P_last_guess * Cr / sqrt(R*T);
    
    [P_res_calc, P_dist_partial] = backwardPressureSolver(P_last_guess,Nt-1,mdot,Cr,R,T,mu2,gamma);
    
    P_dist = [P_dist_partial, P_sump];
    err = P_res_calc - P_reservoir;
end
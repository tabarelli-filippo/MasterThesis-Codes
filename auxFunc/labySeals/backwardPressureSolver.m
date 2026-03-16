function [P_res_calc, P_distribution] = backwardPressureSolver(P_last_cavity,Nt,mdot,Cr,R,T,mu2,gamma)
% BACKPRESSURESOLVER: auxiliary function for labySeals.m
% evaluates pressure distribution backwards, given a guess over pressure in
% the last cavity.

    P_dist_internal = zeros(1, Nt); 
    P_dist_internal(Nt) = P_last_cavity; % Pressure needed before sump
    
    current_P_downstream = P_last_cavity;
    
    % backwards loop over cavities
    for i = (Nt-1):-1:1
        P_upstream = pressure_backward(current_P_downstream, Cr, R, T, mu2, gamma, mdot);
        P_dist_internal(i) = P_upstream;
        current_P_downstream = P_upstream;
    end
    
    P_res_calc = pressure_backward(current_P_downstream, Cr, R, T, mu2, gamma, mdot); 
    P_distribution = [P_res_calc, P_dist_internal]; 
end

function [Pin] = pressure_backward(P_downstream,Cr, R, T, mu2, gamma, mdot)
    Pin_res = @(P_upstream) Pin_wrapperError(P_upstream,P_downstream,Cr,R,T,mu2,gamma,mdot);
    % solution search
    Pin = fzero(Pin_res, [P_downstream*1.01, P_downstream*100]);

end


function res = Pin_wrapperError(Pin, Pout, Cr, R, T, mu2, gamma, mdot)    
    if Pin <= Pout; res = 1e20; return; end

    s = (Pin/Pout)^((gamma-1)/gamma) - 1;
    mu1 = pi / (pi + 2 - 5*s + 2*s^2); % Chaplygin formula

    kin_term = (real(mdot)/mu1/mu2/Cr)^2 * R * T;
    delta_Pressure = Pin^2 - Pout^2;
    
    res = kin_term - delta_Pressure; % residual to minimize
end
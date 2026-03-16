function [P_end_calc, P_distribution] = forwardPressureSolver(P_start, Nt, mdot, Cr, R, T, mu2, gamma, P_targets)
%FORWARDPRESSURESOLVER: auxiliary function for labySeals.m
% evaluates pressure distribution backwards, given a guess over pressure in
% the first cavity 
    P_distribution = zeros(1, Nt + 1); 
    P_distribution(1) = P_start; 
    
    P_upstream = P_start;
    
    % Forward Loop
    for i = 1:Nt
        
        P_guess_local = P_targets(i+1);
        
        P_downstream = pressure_forward(P_upstream, Cr, R, T, mu2, gamma, mdot, P_guess_local);
        
        P_distribution(i+1) = P_downstream;
        P_upstream = P_downstream;
    end
    
    P_end_calc = P_distribution(end);
end

function [Pout] = pressure_forward(P_upstream, Cr, R, T, mu2, gamma, mdot, P_guess_local)
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
        Pout = lb; 
    end
end

function res = Pin_wrapperError_Fwd(Pin, Pout, Cr, R, T, mu2, gamma, mdot)
    
    if Pout >= Pin; res = 1e20; return; end 
    if Pout <= 0; res = -1e20; return; end
    
    s = (Pin/Pout)^((gamma-1)/gamma) - 1;
    mu1 = pi / (pi + 2 - 5*s + 2*s^2);
    
    kin_term = (real(mdot)/mu1/mu2/Cr)^2 * R * T;
    delta_Pressure = Pin^2 - Pout^2;
    
    res = kin_term - delta_Pressure;
end
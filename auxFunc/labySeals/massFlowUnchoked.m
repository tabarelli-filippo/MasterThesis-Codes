function mdot = massFlowUnchoked(P_up, P_down, Cr, R, T, gamma, mu2)
% MASSFLOWUNCHOCKED: auxiliary function for labySeals.m
% evaluates massflow when cavity is not choked
    if P_up <= P_down; mdot = 0; return; end
    
    pr_val = P_down / P_up;
    s = (1/pr_val)^((gamma-1)/gamma) - 1;
    mu1 = pi / (pi + 2 - 5*s + 2*s^2);
    
    mdot = Cr*mu1*mu2*sqrt((P_up^2 - P_down^2)/R/T);
end
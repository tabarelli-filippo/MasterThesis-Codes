function [eta] = viscosityTemp(etaR1,etaR2,tR1,tR2,T)
%VISCOSITYTEMP evaluates viscosity from temperature. Formulation references
% from ASTM D 341 – 03
% INPUT:    etaR1  Viscosity at reference temperature tR1 [Poise]
%           etaR2  Viscosity at reference temperature tR2 [Poise]
%           tR1    Reference Temperature 1 [K]
%           tR1    Reference Temperature 1 [K]
%           T      Evaulation temperature
arguments (Input)
    etaR1 (1,1) double
    etaR2 (1,1) double
    tR1 (1,1) double
    tR2 (1,1) double
    T (1,1) double
end

% coefficient A and B evaluation
alpha = log10(log10([etaR1+ 0.7;etaR2+ 0.7]));
beta = [1, -log10(tR1); 1 -log10(tR2)];

coeff = beta\alpha;
A = coeff(1); B = coeff(2);
% viscosity evaluation
eta = -0.7 + 10.^(10.^(A-B*log10(T)));
end
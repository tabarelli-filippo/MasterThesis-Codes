function [eta] = viscosityTemp(etaR1, etaR2, tR1, tR2, T)
% VISCOSITYTEMP  Estimates the dynamic viscosity of a lubricating oil at a
%   given temperature using the ASTM D 341 viscosity-temperature equation.
%
%   The ASTM D 341 model (also known as the Walther equation) expresses the
%   relationship between kinematic viscosity and absolute temperature on a
%   double-logarithmic scale:
%       log10(log10(nu + 0.7)) = A - B * log10(T)
%
%   where nu is the kinematic viscosity [cSt] and T is the absolute
%   temperature [K]. The constants A and B are determined from two
%   reference viscosity-temperature pairs.
%
%   This function uses dynamic viscosity [Poise = P] as input/output
%   (approximately equal to kinematic viscosity in cSt for oils with
%   density ≈ 1 g/cm³), consistent with the ASTM D 341-03 standard.
%
% SYNTAX
%   eta = viscosityTemp(etaR1, etaR2, tR1, tR2, T)
%
% INPUT ARGUMENTS
%   etaR1 - (scalar double) Dynamic viscosity at reference temperature
%            tR1 [Poise]. Typically the viscosity at 40°C (313.15 K).
%   etaR2 - (scalar double) Dynamic viscosity at reference temperature
%            tR2 [Poise]. Typically the viscosity at 100°C (373.15 K).
%   tR1   - (scalar double) First reference temperature [K]
%   tR2   - (scalar double) Second reference temperature [K]
%   T     - (scalar double) Evaluation temperature [K]
%
% OUTPUT ARGUMENTS
%   eta   - (scalar double) Estimated dynamic viscosity at temperature T
%           [Poise]. Convert to SI units [N·s/m² = Pa·s] by dividing
%           by 10: eta_SI = eta / 10.
%
% THEORY
%   The ASTM D 341 coefficients A and B are found by solving:
%       [1, -log10(tR1)] * [A; B] = log10(log10(etaR1 + 0.7))
%       [1, -log10(tR2)] * [A; B] = log10(log10(etaR2 + 0.7))
%   Then: eta(T) = -0.7 + 10^(10^(A - B*log10(T)))
%
% REFERENCES
%   ASTM D 341-03 (2003). Standard Practice for Viscosity-Temperature
%     Charts for Liquid Petroleum Products. ASTM International.
%
% NOTES
%   - Temperatures must be in Kelvin. Convert from Celsius: T[K] = T[°C] + 273.15.
%   - The ASTM equation is valid for most mineral oils in the range
%     roughly 0°C to 200°C. Extrapolation beyond the reference range
%     should be used with caution.
%   - For use with fluid-film bearing models (bearingMatrix type 7 or
%     nonLinBearingMatrix type 7.1), convert to SI: eta_SI = eta / 10.
%
% EXAMPLE
%   % ISO VG 46 oil: nu = 46 cSt at 40°C, nu = 6.8 cSt at 100°C
%   eta_at_60C = viscosityTemp(46, 6.8, 313.15, 373.15, 333.15);
%   eta_SI = eta_at_60C / 10;   % convert Poise to Pa·s
%   fprintf('Viscosity at 60°C: %.4f Pa·s\n', eta_SI);
%
% SEE ALSO
%   bearingMatrix, nonLinBearingMatrix

arguments (Input)
    etaR1 (1,1) double
    etaR2 (1,1) double
    tR1   (1,1) double
    tR2   (1,1) double
    T     (1,1) double
end

% solve for ASTM D 341 coefficients A and B
alpha = log10(log10([etaR1 + 0.7; etaR2 + 0.7]));
beta  = [1, -log10(tR1); 1, -log10(tR2)];

coeff = beta \ alpha;
A = coeff(1);
B = coeff(2);

% evaluate viscosity at temperature T
eta = -0.7 + 10.^(10.^(A - B*log10(T)));
end

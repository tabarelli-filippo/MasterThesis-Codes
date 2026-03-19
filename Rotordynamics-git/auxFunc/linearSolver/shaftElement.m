function [M0e,C1e,K0e,K1e,beta] = shaftElement(Shaft,L)
% SHAFTELEMENT  Computes the finite element matrices for a uniform circular
%   shaft element used in rotor dynamics analysis.
%
%   Three formulations are available, selected via Shaft.type:
%     Type 1 - Euler-Bernoulli beam: classical slender-beam theory;
%              shear deformation and rotary inertia neglected.
%     Type 2 - Timoshenko beam (Cowper shear coefficient): shear
%              deformation and rotary inertia included; recommended for
%              short or stubby shaft sections.
%     Type 3 - Timoshenko beam (Hutchinson shear coefficient): as Type 2
%              but uses the Hutchinson (2001) formula for the shear
%              correction factor, which is more accurate for hollow shafts.
%
%   Optional effects that can be superimposed:
%     - Axial force (geometric stiffness, tension/compression)
%     - Torque (follower-force geometric stiffness)
%     - Internal (material/hysteretic) damping via the beta coefficient
%
%   The element has 8 DOFs: [u1, v1, theta_u1, theta_v1, u2, v2, theta_u2, theta_v2]
%   where 1 and 2 denote the two end nodes.
%
% SYNTAX
%   [M0e, C1e, K0e, K1e, beta] = shaftElement(Shaft, L)
%
% INPUT ARGUMENTS
%   Shaft - (struct) Shaft element definition. Required fields:
%             .type    - element type (1, 2, or 3)
%             .d_ext1  - external diameter [m]
%             .d_int1  - internal diameter [m] (0 for solid shaft)
%             .rho     - material density [kg/m³]
%             .E       - Young's modulus [Pa]
%             .G       - shear modulus [Pa]
%           Optional fields:
%             .AxialForce - axial force along the element [N]
%                           (positive = tension, negative = compression)
%             .Torque     - torque applied along the element [N·m]
%             .beta       - internal damping coefficient [s]
%                           (proportional damping ratio: C = beta * K)
%   L     - (double, scalar) Element length [m]
%
% OUTPUT ARGUMENTS
%   M0e  - (8x8 double) Element consistent mass matrix [kg, kg·m²].
%          Includes translational and (if type 2 or 3) rotary inertia.
%   C1e  - (8x8 double) Element gyroscopic matrix [kg·m²].
%          Assembled as Omega * C1e in the global EOM.
%   K0e  - (8x8 double) Element elastic stiffness matrix [N/m, N·m/rad].
%          Includes shear, axial-force, and torque effects if specified.
%   K1e  - (8x8 double) Speed-dependent stiffness matrix due to internal
%          damping [N/m]. Non-zero only when beta > 0. Assembled as
%          Omega * K1e; skew-symmetric, hence destabilizing.
%   beta - (double) Internal damping coefficient [s] as read from the
%          Shaft struct (0 if not specified).
%
% THEORY
%   For a hollow circular cross-section:
%     A  = pi/4 * (d_ext^2 - d_int^2)    cross-sectional area
%     J  = pi/64 * (d_ext^4 - d_int^4)   second moment of area
%   Shear correction factor kappa (Cowper 1966 or Hutchinson 2001) is
%   used to compute the shear coefficient phi = 12*E*J/(G*kappa*A*L^2).
%
% REFERENCES
%   Cowper, G.R. (1966). The shear coefficient in Timoshenko's beam theory.
%     Journal of Applied Mechanics, 33(2), 335-340.
%   Hutchinson, J.R. (2001). Shear coefficients for Timoshenko beam theory.
%     Journal of Applied Mechanics, 68(1), 87-92.
%
% EXAMPLE
%   Shaft.type  = 2;       % Timoshenko
%   Shaft.d_ext1 = 0.08;   % [m]
%   Shaft.d_int1 = 0.0;    % solid
%   Shaft.rho   = 7800;    % [kg/m^3]
%   Shaft.E     = 2.1e11;  % [Pa]
%   Shaft.G     = 8.1e10;  % [Pa]
%   Shaft.beta  = 1e-5;    % internal damping [s]
%   L = 0.25;              % [m]
%   [M0e, C1e, K0e, K1e, beta] = shaftElement(Shaft, L);
%
% SEE ALSO
%   shaftTaperEl, diskElement, rotorMatrix

if isfield(Shaft, 'AxialForce')
    AxialForce = Shaft.AxialForce;
else
    AxialForce = 0;
end

if isfield(Shaft, 'Torque')
    Torque = Shaft.Torque;
else
    Torque = 0;
end

if isfield(Shaft, 'beta')
    beta = Shaft.beta;
else
    beta = 0;
end

Shaft_Type = Shaft.type;
d_ext = Shaft.d_ext1;
d_int = Shaft.d_int1;
rho = Shaft.rho; 
E = Shaft.E;
G = Shaft.G;

J = pi / 64 * (d_ext^4 - d_int^4);
A = pi / 4 *(d_ext^2-d_int^2);


switch Shaft_Type
    case 1 % Euler-Bernoulli Beam
        include_shear_effects = 0;
        include_rotary_inertia = 0;
        include_gyroscopic = 1;
    case {2,3} % Timoshenko Beam
        include_shear_effects = 1;
        include_rotary_inertia = 1;
        include_gyroscopic = 1;
end

% Including Shear Effect
if (include_shear_effects~=0)
   Poisson = 0.5*(E/G) - 1;
   r = d_int/d_ext;
   r2 = r^2; 
   r12 = (1+r2)^2;
   if Shaft_Type == 3
        kappa_num = 6 * r12 * (1+Poisson)^2;
        kappa_den = r12*(7+12*Poisson+4*Poisson^2) + 4*r2*(5+6*Poisson+2*Poisson^2);
        Kappa = kappa_num/kappa_den;  %Hutchinson
   else
        Kappa = 6*r12*(1+Poisson)/(r12*(7+6*Poisson)+r2*(20+12*Poisson));  %Cowper
   end 
   shear_coeff = 12*E*J/(G*Kappa*A*L*L);
else
   shear_coeff = 0.0;
end


phi = shear_coeff;

%% Stiffness Element Matrix
K0e = [12     0           0         6*L  -12     0           0         6*L;
        0    12        -6*L           0    0   -12        -6*L           0;
        0  -6*L (4+phi)*L*L           0    0   6*L (2-phi)*L*L           0;
      6*L     0           0 (4+phi)*L*L -6*L     0           0 (2-phi)*L*L;
      -12     0           0        -6*L   12     0           0        -6*L;
        0   -12         6*L           0    0    12         6*L           0;
        0  -6*L (2-phi)*L*L           0    0   6*L (4+phi)*L*L           0;
      6*L     0           0 (2-phi)*L*L -6*L     0           0 (4+phi)*L*L];
K0e = E*J*K0e/( (1+phi)*L^3 );

%% Mass Element Matrix
m1 = 312 + 588*phi + 280*phi^2;
m2 = (44 + 77*phi + 35*phi^2)*L;
m3 = 108 + 252*phi + 140*phi^2;
m4 = -(26 + 63*phi + 35*phi^2)*L;
m5 = (8 + 14*phi +7*phi*phi)*L^2;
m6 = -(6 + 14*phi +7*phi^2)*L^2;

M0e = [ m1    0    0   m2   m3    0    0   m4;
         0   m1  -m2    0    0   m3  -m4    0;
         0  -m2   m5    0    0   m4   m6    0;
        m2    0    0   m5  -m4    0    0   m6;
        m3    0    0  -m4   m1    0    0  -m2; 
         0   m3   m4    0    0   m1   m2    0;
         0  -m4   m6    0    0   m2   m5    0;
        m4    0    0   m6  -m2    0    0   m5];
M0e = rho*A*L*M0e/(840*(1+phi)^2);

% include the rotary inertia effects in the mass matrix
if (include_rotary_inertia~=0)
   phi = shear_coeff;
   if (include_shear_effects==0)
       phi = 0; 
   end

   m7 = 36;
   m8 = (3 - 15*phi)*L;
   m9 = (4 + 5*phi +10*phi^2)*L^2;
   m10 = (-1 - 5*phi + 5*phi^2)*L^2;
   Ms =  [ m7    0    0   m8  -m7    0    0   m8;
            0   m7  -m8    0    0  -m7  -m8    0;
            0  -m8   m9    0    0   m8  m10    0;
           m8    0    0   m9  -m8    0    0  m10;
          -m7    0    0  -m8   m7    0    0  -m8;
            0  -m7   m8    0    0   m7   m8    0;
            0  -m8  m10    0    0   m8   m9    0;
           m8    0    0  m10  -m8    0    0   m9];
   Ms = rho*J*Ms/(30*L*(1+phi)^2);
   M0e = M0e + Ms;
end

%%  Gyroscopic Element Matrix
if (include_gyroscopic==0)
   C1e = zeros(8,8);
else
   phi = shear_coeff;
   if (include_shear_effects==0)
       phi = 0; 
   end
   g1 = 36;
   g2 = (3-15*phi)*L;
   g3 = (4 + 5*phi + 10*phi^2)*L^2;
   g4 = (-1 - 5*phi + 5*phi^2)*L^2;
   C1e = [  0  -g1   g2    0    0   g1   g2    0;
           g1    0    0   g2  -g1    0    0   g2;   
          -g2    0    0  -g3   g2    0    0  -g4;
            0  -g2   g3    0    0   g2   g4    0;
            0   g1  -g2    0    0  -g1  -g2    0;
          -g1    0    0  -g2   g1    0    0  -g2;
          -g2    0    0  -g4   g2    0    0  -g3;
            0  -g2   g4    0    0   g2   g3    0];
   C1e = - rho*J*C1e/(15*L*(1+phi)^2);
end
%%  Axial force stiffness Matrix
if AxialForce ~= 0
   phi = shear_coeff; 
   if (include_shear_effects==0)
       phi = 0;
   end

   k1 = 72 + 120*phi + 60*phi^2;
   k2 = 6*L;
   k3 = (8 + 10*phi + 5*phi^2)*L^2;
   k4 = (-2 - 10*phi - 5*phi^2)*L^2;
   KFe = [ k1    0    0   k2  -k1    0    0   k2;
            0   k1  -k2    0    0  -k1  -k2    0;
            0  -k2   k3    0    0   k2   k4    0;
           k2    0    0   k3  -k2    0    0   k4;
          -k1    0    0  -k2   k1    0    0  -k2; 
            0  -k1   k2    0    0   k1   k2    0;
            0  -k2   k4    0    0   k2   k3    0;
           k2    0    0   k4  -k2    0    0   k3];
   KFe = AxialForce*KFe/(60*L*(1+phi)^2);
   K0e = K0e + KFe;

end
%%  Torque stiffness Matrix
% this effect is based on Euler-Bernoulli formulation
if Torque ~= 0
   KTe = [  0    0    1    0    0    0   -1    0;
            0    0    0    1    0    0    0   -1;
            1    0    0 -L/2   -1    0    0  L/2;
            0    1  L/2    0    0   -1 -L/2    0;
            0    0   -1    0    0    0    1    0; 
            0    0    0   -1    0    0    0    1;
           -1    0    0 -L/2    1    0    0  L/2;
            0   -1  L/2    0    0    1 -L/2    0];
   KTe = Torque*KTe/L;
   K0e = K0e + KTe;
end

%% Internal Damping / Material Damping
if beta ~= 0
    K1e = [ 0    12         -6*L            0    0   -12         -6*L            0;
        -12     0            0         -6*L   12     0            0         -6*L;
        6*L     0            0  (4+phi)*L*L -6*L     0            0  (2-phi)*L*L;
        0   6*L -(4+phi)*L*L            0    0  -6*L -(2-phi)*L*L            0;
        0   -12          6*L            0    0    12          6*L            0;
        12     0            0          6*L  -12     0            0          6*L;
        6*L     0            0  (2-phi)*L*L -6*L     0            0  (4+phi)*L*L;
        0   6*L -(2-phi)*L*L            0    0  -6*L -(4+phi)*L*L            0];
    K1e = E*J*K1e/( (1+phi)*L^3 );
else
    K1e = zeros(8,8);
end
end

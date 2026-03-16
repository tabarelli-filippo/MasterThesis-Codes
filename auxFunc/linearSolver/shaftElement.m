function [M0e,C1e,K0e,K1e,beta] = shaftElement(Shaft,L)
%SHAFTELEMENT computes the element matrix for a given shaft description.
% Only circular shaft. For each element see Manual
%
%INPUT: Shaft   Structure
%       L       Element length
%       optional: 
%       AxialForce  Axial force along the element
%       Torque      Torque effect along the element
%       beta        Internal Damping / Material Damping 
%
%OUTPUT:K0e  stiffness matrix
%       C1e  gyroscopic matrix
%       M0e  mass matrix
%       K1e  speed dependent contribution to the  stiffness matrix due to the 
%            internal damping

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
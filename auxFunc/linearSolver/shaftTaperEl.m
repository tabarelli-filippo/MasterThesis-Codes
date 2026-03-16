function [M0e,C1e,K0e,K1e] = shaftTaperEl(Shaft,L)
%SHAFTELEMENT computes the element matrix for a given shaft description.
% Only tapered shaft. For each element see Manual
%
%INPUT: Shaft   Structure
%       L       Element length
%
%OUTPUT:K0e  stiffness matrix
%       C1e  gyroscopic matrix
%       M0e  mass matrix
%       K1e  speed dependent contribution to the  stiffness matrix due to the 
%            internal damping

Shaft_Type = Shaft.type - 20;

if isfield(Shaft, 'AxialForce')
    AxialForce = Shaft.AxialForce;
else
    AxialForce = 0;
end

%node 1
doj = Shaft.d_ext1;
dij = Shaft.d_int1;
%node 2
dok = Shaft.d_ext2;
dik = Shaft.d_int2;

rho = Shaft.rho; 
E = Shaft.E;
G = Shaft.G;


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

% define intermediate constants
nu = 0.5*(E/G) - 1;
rok = dok/2; rik = dik/2; roj = doj/2; rij = dij/2;
romean = (rok+roj)/2; rimean = (rik+rij)/2; 
mu2 = (rimean/romean)^2;

Am = pi*(romean^2-rimean^2); 
Im = pi*(romean^4-rimean^4)/4;

if include_shear_effects~=0
    if Shaft_Type == 3
        kappa_num = 6 * (1+mu2)^2 * (1+nu)^2;
        kappa_den = (1+mu2)^2*(7+12*nu+4*nu^2) + 4*mu2*(5+6*nu+2*nu^2);
        kappa = kappa_num/kappa_den;  %Hutchinson 2001
    else 
        kappa = 6*(1+nu)*(1+mu2)^2/((7+6*nu)*(1+mu2)^2+(20+12*nu)*mu2); %Cowper 1966
    end
    chi = 1/kappa;
    phi = 12*E*Im*chi/(G*Am*L^2);
else 
    phi = 0;
end

Dro = rok-roj; Dri = rik-rij;
Aj = pi*(doj^2-dij^2)/4; 
Itj = pi*(doj^4-dij^4)/64;


a1 = 2*pi*(roj*Dro-rij*Dri)/Aj;
b1 = pi*(Dro^2-Dri^2)/Aj;
a2 = pi*(roj^3*Dro-rij^3*Dri)/Itj;
b2 = 3*pi*(roj^2*Dro^2-rij^2*Dri^2)/(2*Itj);
g2 = pi*(roj*Dro^3-rij*Dri^3)/Itj;
d2 = pi*(Dro^4-Dri^4)/(4*Itj);


%% Stiffness Element Matrix
k1 = 1260+630*a2+504*b2+441*g2+396*d2;
k2 = L*(630+210*a2+147*b2+126*g2+114*d2-phi*(105*a2+105*b2+94.5*g2+84*d2));
k3 = L*(630+420*a2+357*b2+315*g2+282*d2+phi*(105*a2+105*b2+94.5*g2+84*d2));
k4 = L^2*(420+210*phi+105*phi^2+a2*(105+52.5*phi^2)+b2*(56-35*phi+35*phi^2) ...
      +g2*(42-42*phi+26.25*phi^2)+d2*(36-42*phi+21*phi^2));
k5 = L^2*(210-210*phi-105*phi^2+a2*(105-105*phi-52.5*phi^2)+b2*(91-70*phi-35*phi^2) ...
      +g2*(84-52.5*phi-26.25*phi^2)+d2*(78-42*phi-21*phi^2));
k6 = L^2*(420+210*phi+105*phi^2+a2*(315+210*phi+52.5*phi^2)+b2*(266+175*phi+35*phi^2) ...
      +g2*(231+147*phi+26.25*phi^2)+d2*(204+126*phi+21*phi^2));
k7 = 12+6*a1+4*b1;
k8 = L*(6+3*a1+2*b1);
k9 = L^2*(3+1.5*a1+b1);

ka = [k1   0   0  k2 -k1   0   0  k3;
       0  k1 -k2   0   0 -k1 -k3   0;
       0 -k2  k4   0   0  k2  k5   0;
      k2   0   0  k4 -k2   0   0  k5;
     -k1   0   0 -k2  k1   0   0 -k3;
       0 -k1  k2   0   0  k1  k3   0;
       0 -k3  k5   0   0  k3  k6   0;
      k3   0   0  k5 -k3   0   0  k6];

kb = [k7   0   0  k8 -k7   0   0  k8;
       0  k7 -k8   0   0 -k7 -k8   0;
       0 -k8  k9   0   0  k8  k9   0;
      k8   0   0  k9 -k8   0   0  k9;
     -k7   0   0 -k8  k7   0   0 -k8;
       0 -k7  k8   0   0  k7  k8   0;
       0 -k8  k9   0   0  k8  k9   0;
      k8   0  0   k9 -k8   0   0  k9];

ke = E*Itj/(105*L^3*(1+phi)^2)*(ka+105*phi*kb);
K0e = ke;
K1e = zeros(8,8); 
%% axial force
if AxialForce ~= 0
  k10 = 36+60*phi+3*phi^2;
  k11 = L*3;
  k12 = L^2*(4+5*phi+2.5*phi^2);
  k13 = L^2*(1+5*phi+2.5*phi^2);

  kG = AxialForce/(30*L*(1+phi)^2)*...
     [k10    0    0  k11 -k10    0    0  k11;
        0  k10 -k11    0    0 -k10 -k11    0;
        0 -k11  k12    0    0  k11 -k13    0; 
      k11    0    0  k12 -k11    0    0 -k13; 
     -k10    0    0 -k11  k10    0    0 -k11;
        0 -k10  k11    0    0  k10  k11    0;
        0 -k11 -k13    0    0  k11  k12    0;
      k11    0    0 -k13 -k11    0    0  k12];    
  K0e = K0e + kG;
end


%% Mass Element Matrix
% Mass matrix
m1 = (468+882*phi+420*phi^2)+a1*(108+210*phi+105*phi^2)+b1*(38+78*phi+42*phi^2);
m2 = L*((66+115.5*phi+52.5*phi^2)+a1*(21+40.5*phi+21*phi^2)+b1*(8.5+18*phi+10.5*phi^2));
m3 = (162+378*phi+210*phi^2)+a1*(81+189*phi+105*phi^2)+b1*(46+111*phi+63*phi^2);
m4 = L*((39+94.5*phi+52.5*phi^2)+a1*(18+40.5*phi+21*phi^2)+b1*(9.5+21*phi+10.5*phi^2));
m5 = L^2*((12+21*phi+10.5*phi^2)+a1*(4.5+9*phi+5.25*phi^2)+b1*(2+4.5*phi+3*phi^2));
m6 = L*((39+94.5*phi+52.5*phi^2)+a1*(21+54*phi+31.5*phi^2)+b1*(12.5+34.5*phi+21*phi^2));
m7 = L^2*((9+21*phi+10.5*phi^2)+a1*(4.5+10.5*phi+5.25*phi^2)+b1*(2.5+6*phi+3*phi^2));
m8 = (468+882*phi+420*phi^2)+a1*(360+672*phi+315*phi^2)+b1*(290+540*phi+252*phi^2);
m9 = L*((66+115.5*phi+52.5*phi^2)+a1*(45+75*phi+31.5*phi^2)+b1*(32.5+52.5*phi+21*phi^2));
m10 = L^2*((12+21*phi+10.5*phi^2)+a1*(7.5+12*phi+5.25*phi^2)+b1*(5+7.5*phi+3*phi^2));

mT = rho*Aj*L/(1260*(1+phi)^2)* ...
   [m1   0   0  m2  m3   0   0 -m4;
     0  m1 -m2   0   0  m3  m4   0;
     0 -m2  m5   0   0 -m6 -m7   0;
    m2   0   0  m5  m6   0   0 -m7;
    m3   0   0  m6  m8   0   0 -m9;
     0  m3 -m6   0   0  m8  m9   0;
     0  m4 -m7   0   0  m9 m10   0;
   -m4   0   0 -m7 -m9   0   0 m10]; 
M0e = mT;

m11 = 252+126*a2+72*b2+45*g2+30*d2;
m12 = L*(21-105*phi+a2*(21-42*phi)+b2*(15-21*phi)+g2*(10.5-12*phi)+d2*(7.5-7.5*phi));
m13 = L*(21-105*phi+a2*(-63*phi)-b2*(6+42*phi)-g2*(7.5+30*phi)-d2*(7.5+22.5*phi));
m14 = L^2*((28+35*phi+70*phi^2)+a2*(7-7*phi+17.5*phi^2)+b2*(4-7*phi+7*phi^2) ...
       +g2*(2.75-5*phi+3.5*phi^2)+d2*(2-3.5*phi+2*phi^2));
m15 = L^2*((7+35*phi-35*phi^2)+a2*(3.5+17.5*phi-17.5*phi^2)+b2*(3+10.5*phi-10.5*phi^2) ...
       +g2*(2.75+7*phi-7*phi^2)+d2*(2.5+5*phi-5*phi^2));
m16 = L^2*((28+35*phi+70*phi^2)+a2*(21+42*phi+52.5*phi^2)+b2*(18+42*phi+42*phi^2) ...
       +g2*(16.25+40*phi+35*phi^2)+d2*(15+37.5*phi+30*phi^2));

%% Include rotary inertia
if include_rotary_inertia ~=0
  mR = rho*Itj/(210*L*(1+phi)^2)*...
     [m11    0    0  m12 -m11    0    0  m13;
        0  m11 -m12    0    0 -m11 -m13    0;
        0 -m12  m14    0    0  m12 -m15    0;
      m12    0    0  m14 -m12    0    0 -m15; 
     -m11    0    0 -m12  m11    0    0 -m13;
        0 -m11  m12    0    0  m11  m13    0;
        0 -m13 -m15    0    0  m13  m16    0;
      m13    0    0 -m15 -m13    0    0  m16];      
  M0e = mT + mR;
end 
%%  Gyroscopic Element Matrix
if include_gyroscopic ~=0
    ge = rho*Itj/(105*L*(1+phi)^2)*...
        [  0 -m11  m12    0    0  m11  m13    0
        m11    0    0  m12 -m11    0    0  m13
        -m12    0    0 -m14  m12    0    0  m15
        0 -m12  m14    0    0  m12 -m15    0
        0  m11 -m12    0    0 -m11 -m13    0
        -m11    0    0 -m12  m11    0    0 -m13
        -m13    0    0  m15  m13    0    0 -m16
        0 -m13 -m15    0    0  m13  m16    0];
    C1e = - ge;
else
    C1e = zeros(8,8);
end

end
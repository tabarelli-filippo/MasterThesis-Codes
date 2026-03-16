function [time,q,qdot,speed] = runUp(Rotor,t_span,alpha,initPos,initVel,options)
%TIMESIMULATION computes time response solving non-linear ODE associated
%with the problem. Might require some time
%
%INPUT: Rotor       structure
%       t_span      t_span = [t_init,t_end] initial and endind time of the 
%                   simulation. [s]
%       initPos     initial condition on position
%       initVel     initial condition on velocity
%       alpha       rotor angle time law: 
%                   [a2 a1 a0] phi = a2*t^2 + a1*t + a0
%
%OPTIONAL:  nr      number of degrees reduced, if nr = 0, no reduction is
%                   applied (default). Model reduction is based on undamped
%                   modes
arguments (Input)
    Rotor (1,1) struct
    t_span (1,2) double
    alpha (:,:) double {mustBeVector}
    initPos (:,1) double
    initVel (:,1) double
    options.nr (1,1) double = 0;
end
nr = options.nr;
n_nodes = numel(Rotor.nodes);
ndof = 4 * n_nodes;

% active DoFs
rotSpeedref = 0.0;
[~,~,~,zero_dof,~] = bearingMatrix(Rotor,rotSpeedref);
dof = 1:ndof;
dof(zero_dof) = [];
nzero = length(zero_dof);
ncdof = ndof - nzero;

% rotor matrices
[M0,C0,C1,K0,K1] = rotorMatrix(Rotor);
% bearing static matrices
rotSpeedref = 0.0; 
[Mb_st, Cb_st, Kb_st, zero_dof, ~] = bearingMatrix(Rotor, rotSpeedref);

M_tot = M0 + Mb_st;
K_tot = K0 + Kb_st;
C_tot = C0 + Cb_st;
%% Model reduction
if nr > 0 && nr < ncdof
   [T, ~] = reductionMatrix(Rotor,rotSpeedref,nr);

   Mr = T.' * M_tot(dof,dof) * T;
   Kr = T.' * K_tot(dof,dof) * T;
   Cr = T.' * C_tot(dof,dof) * T;
   
   C1r = T.' * C1(dof,dof) * T; 
   K1r = T.' * K1(dof,dof) * T;
   
   q0 = T.' * initPos(dof);
   qdot0 = T.' * initVel(dof);
else
   nr = ncdof;
   T = eye(ncdof);
   
   Mr = M_tot(dof,dof);
   Kr = K_tot(dof,dof);
   Cr = C_tot(dof,dof);
   C1r = C1(dof,dof);
   K1r = K1(dof,dof);
   
   q0 = initPos(dof);
   qdot0 = initVel(dof);
end
% static A matrix component pre-evalutation
A21_st = -Mr\Kr;
A22_st = -Mr\Cr;

% gyro A matrix component pre-evalutation
A21_dyn = -Mr\K1r;
A22_dyn= -Mr\C1r;

% static forcing
T_force = Mr\T.';
[unForce, bendForce] = forcing(Rotor);
 
initCond = [q0; qdot0];
%% Equation definition
ode_fun = @(t, y) system(t, y, A21_st, A22_st, A21_dyn, A22_dyn, unForce,...
    bendForce, T_force, Rotor, alpha, T, dof, nr, ndof);

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

[time,y] = ode15s(ode_fun, t_span, initCond, options);
% speed over time
speed = 2*alpha(1)*time + alpha(2)*ones(size(time));

q = zeros(length(time), ndof);
qdot = zeros(length(time), ndof);

if nr < ncdof
    q(:, dof) = y(:, 1:nr) * T.';
    qdot(:, dof) = y(:, nr+1:end) * T.';
else
    q(:, dof) = y(:, 1:nr);
    qdot(:, dof) = y(:, nr+1:end);
end
end

%% ODE FUNCTION
function [ydot] = system(t, y, A21_st, A22_st, A21_dyn, A22_dyn, unForce,...
    bendForce, T_force, Rotor, alpha, T, dof, nr, ndof)

    % runup
    phi = alpha(1)*t*t + alpha(2)*t + alpha(3);
    phidot = 2*alpha(1)*t + alpha(2);
    phidotdot =  2*alpha(1);
    
    A21 = A21_st + phidot*A21_dyn;
    A22 = A22_st + phidot*A22_dyn;

    % position and velocity estraction
    q_red = y(1:nr);
    q_dot_red = y(nr+1:end);
    
    q = zeros(ndof, 1);
    q_dot = zeros(ndof, 1);
    q(dof) = T * q_red;
    q_dot(dof) = T * q_dot_red;
    
    % nonLinear Forcing
    [nonLinBearForce] = nonLinBearingMatrix(Rotor,phidot,q,q_dot);
    
    % forcing
    rot_factor = exp(1i * phi); % e^(j*Omega*t)
    ff = bendForce + unForce*(phidot*phidot - 1i*phidotdot); 
    Force = real(ff*rot_factor) + nonLinBearForce;
    
    ForceR = T_force*Force(dof);

    % variable state matrices definition
    A = [zeros(nr,nr) eye(nr); A21, A22];
    F = [zeros(nr,1);    ForceR];
    
    ydot = A*y + F;
end
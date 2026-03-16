function [time,q,qdot] = timeSimulation(Rotor,t_span,rotorSpeed,initPos,initVel,options)
%TIMESIMULATION computes time response solving non-linear ODE associated
%with the problem. Might require some time
%
%INPUT: Rotor       structure
%       t_span      t_span = [t_init,t_end] initial and endind time of the 
%                   simulation. [s]
%       rotorSpeed  Speed evaluated = constant. NO RUN-UP
%       initPos     initial condition on position
%       initVel     initial condition on velocity
%OPTIONAL:  nr      number of degrees reduced, if nr = 0, no reduction is
%                   applied (default). Model reduction is based on undamped
%                   modes
arguments (Input)
    Rotor (1,1) struct
    t_span (1,2) double
    rotorSpeed (1,1) double
    initPos (:,1) double
    initVel (:,1) double
    options.nr (1,1) double = 0;
end
nr = options.nr;
n_nodes = numel(Rotor.nodes);
ndof = 4 * n_nodes;
% active DoFs
[~,~,~,zero_dof,~] = bearingMatrix(Rotor,rotorSpeed);
dof = 1:ndof;
dof(zero_dof) = [];
nzero = length(zero_dof);
ncdof = ndof - nzero;

[M0,C0,C1,K0,K1] = rotorMatrix(Rotor);
[Mb,Cb,Kb,~,~]   = bearingMatrix(Rotor,rotorSpeed);

% global matrices
M = M0 + Mb;
K = K0 + Kb + rotorSpeed*K1;
C = C0 + Cb + rotorSpeed*C1;


%% Model reduction
if nr > 0 && nr < ncdof
   [T, ~] = reductionMatrix(Rotor,rotorSpeed,nr); 
   
   Mr = T.' * M(dof,dof) * T;
   Cr = T.' * C(dof,dof) * T;
   Kr = T.' * K(dof,dof) * T;
   
   q0 = T.' * initPos(dof);
   qdot0 = T.' * initVel(dof);
else
   nr = ncdof;
   T = eye(ncdof); 

   Mr = M(dof,dof);
   Cr = C(dof,dof);
   Kr = K(dof,dof);
   
   q0 = initPos(dof);
   qdot0 = initVel(dof);
end
% variable state matrix definition
A = [zeros(nr,nr),eye(nr); -Mr\Kr, -Mr\Cr];
% initial condition evaluation
initCond = [q0;qdot0];
% static forcing
[unForce, bendForce,gravityAcc] = forcing(Rotor);

gravityForce = M * gravityAcc;
Force_static = bendForce + unForce*rotorSpeed^2;
T_force = Mr\T.';

%% Equation definition
ode_fun = @(t, y) system(t, y, A, T_force, Force_static, gravityForce, Rotor, rotorSpeed, T, dof, nr, ndof);
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);

[time,y] = ode15s(ode_fun, t_span, initCond, options);

if nr < ncdof
    q_red = y(:, 1:nr);
    qdot_red = y(:, nr+1:end);
    
    q = zeros(length(time), ndof);
    qdot = zeros(length(time), ndof);
    
    q(:, dof) = q_red * T.';
    qdot(:, dof) = qdot_red * T.';
else
    q = zeros(length(time), ndof);
    qdot = zeros(length(time), ndof);
    q(:, dof) = y(:, 1:nr);
    qdot(:, dof) = y(:, nr+1:end);
end

end

%% ODE FUNCTION
function [ydot] = system(t, y, A, T_force, Force_static, gravityForce, Rotor, rotorSpeed, T, dof, nr, ndof)      
    % position and velocity extraction for matrices
    % reduced positions
    q_red = y(1:nr);
    q_dot_red = y(nr+1:end);

    q = zeros(ndof, 1);
    q_dot = zeros(ndof, 1);

    q(dof) = T * q_red;
    q_dot(dof) = T * q_dot_red;
   
    % non-linear forces
    [nonLinBearForce] = nonLinBearingMatrix(Rotor,rotorSpeed,q,q_dot);

    % forcing
    rot_factor = exp(1i * rotorSpeed * t); % e^(jOmega*t)
    Force = real(Force_static*rot_factor) + nonLinBearForce + gravityForce;
    ForceR = T_force*Force(dof);
    % variable state matrices definition
    F = [zeros(nr,1);   ForceR];

    ydot = A*y + F;
end
function [time, q, qdot] = timeSimulation(Rotor, t_span, rotorSpeed, initPos, initVel, options)
% TIMESIMULATION  Computes the time-domain nonlinear response of a rotor-
%   bearing system at a constant rotational speed.
%
%   The full nonlinear equations of motion are integrated using ode15s
%   (stiff solver). The right-hand side includes:
%     - Linear speed-dependent system matrices (M, C + Omega*C1, K + Omega*K1)
%     - Nonlinear bearing forces (via nonLinBearingMatrix)
%     - Constant synchronous forcing: unbalance, shaft bow (rotated by
%       e^{i*Omega*t}) and gravity
%
%   The state-space formulation is:
%       [q_dot ]   [ 0    I  ] [q  ]   [       0          ]
%       [q_ddot] = [-M\K -M\C] [q_dot] + [M\(F_lin + F_nonlin)]
%
%   Optional modal reduction (using reductionMatrix) can be applied to
%   reduce computational cost for large models.
%
% SYNTAX
%   [time, q, qdot] = timeSimulation(Rotor, t_span, rotorSpeed, initPos, initVel)
%   [time, q, qdot] = timeSimulation(Rotor, t_span, rotorSpeed, initPos, initVel, "nr", nr)
%
% INPUT ARGUMENTS
%   Rotor      - (1x1 struct) Rotor data structure (see rotorMatrix,
%                bearingMatrix, forcing for field definitions)
%   t_span     - (1x2 double) [t_start, t_end] integration interval [s]
%   rotorSpeed - (scalar double) Constant shaft speed [rad/s]. All
%                speed-dependent matrices are evaluated at this value.
%   initPos    - (ndof x 1 double) Initial nodal displacement vector [m].
%                Use zeros(ndof,1) for a start from the undeflected shape.
%   initVel    - (ndof x 1 double) Initial nodal velocity vector [m/s].
%                Use zeros(ndof,1) for a start from rest.
%
% NAME-VALUE OPTIONS
%   "nr" - (non-negative integer, default 0) Number of modes for modal
%          reduction. nr = 0 runs the full physical model. Choose nr such
%          that the retained bandwidth (from reductionMatrix) covers the
%          highest frequency of interest.
%
% OUTPUT ARGUMENTS
%   time - (ntime x 1 double) Time vector returned by the ODE solver [s]
%   q    - (ntime x ndof double) Nodal displacement history [m].
%          Rows = time steps, columns = DOFs (zero at constrained DOFs).
%   qdot - (ntime x ndof double) Nodal velocity history [m/s].
%
% NOTES
%   - ODE solver tolerances: RelTol = 1e-3, AbsTol = 1e-6 (ode15s).
%   - Gravity is included via a constant force vector M * g_acc.
%   - Steady-state response can be extracted from the tail of q (after
%     transients have decayed), or by using FRF for linear systems.
%   - ndof = 4 * n_nodes
%
% EXAMPLE
%   ndof = 4 * numel(Rotor.nodes);
%   Omega = 1000 * pi/30;    % 1000 rpm in rad/s
%   [t, q, qd] = timeSimulation(Rotor, [0, 2], Omega, ...
%                               zeros(ndof,1), zeros(ndof,1), "nr", 6);
%   % Extract orbit at node 3
%   node = 3;
%   poincareMap(t, q(:, 4*node-3), q(:, 4*node-2), Omega);
%
% SEE ALSO
%   runUp, reductionMatrix, nonLinBearingMatrix, forcing, FRF,
%   poincareMap, plotBifurcation, DFT

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

% active (unconstrained) DOFs
[~,~,~,zero_dof,~] = bearingMatrix(Rotor,rotorSpeed);
dof = 1:ndof;
dof(zero_dof) = [];
nzero = length(zero_dof);
ncdof = ndof - nzero;

% global system matrices (speed-dependent)
[M0,C0,C1,K0,K1] = rotorMatrix(Rotor);
[Mb,Cb,Kb,~,~]   = bearingMatrix(Rotor,rotorSpeed);

M = M0 + Mb;
K = K0 + Kb + rotorSpeed*K1;
C = C0 + Cb + rotorSpeed*C1;

%% Optional modal reduction
if nr > 0 && nr < ncdof
   [T, ~] = reductionMatrix(Rotor, rotorSpeed, nr); 
   
   Mr = T.' * M(dof,dof) * T;
   Cr = T.' * C(dof,dof) * T;
   Kr = T.' * K(dof,dof) * T;
   
   q0    = T.' * initPos(dof);
   qdot0 = T.' * initVel(dof);
else
   nr = ncdof;
   T = eye(ncdof); 

   Mr = M(dof,dof);
   Cr = C(dof,dof);
   Kr = K(dof,dof);
   
   q0    = initPos(dof);
   qdot0 = initVel(dof);
end

% constant state-space matrix
A = [zeros(nr,nr), eye(nr); -Mr\Kr, -Mr\Cr];

initCond = [q0; qdot0];

% forcing
[unForce, bendForce, gravityAcc] = forcing(Rotor);
gravityForce  = M * gravityAcc;
Force_static  = bendForce + unForce*rotorSpeed^2;
T_force = Mr\T.';

%% ODE integration
ode_fun = @(t, y) system(t, y, A, T_force, Force_static, gravityForce, ...
    Rotor, rotorSpeed, T, dof, nr, ndof);
ode_opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);

[time, y] = ode15s(ode_fun, t_span, initCond, ode_opts);

% reconstruct full-DOF physical solution
if nr < ncdof
    q_red    = y(:, 1:nr);
    qdot_red = y(:, nr+1:end);
    
    q    = zeros(length(time), ndof);
    qdot = zeros(length(time), ndof);
    
    q(:, dof)    = q_red * T.';
    qdot(:, dof) = qdot_red * T.';
else
    q    = zeros(length(time), ndof);
    qdot = zeros(length(time), ndof);
    q(:, dof)    = y(:, 1:nr);
    qdot(:, dof) = y(:, nr+1:end);
end

end

%% -------------------------------------------------------------------------
%  ODE RIGHT-HAND SIDE
% -------------------------------------------------------------------------
function [ydot] = system(t, y, A, T_force, Force_static, gravityForce, ...
    Rotor, rotorSpeed, T, dof, nr, ndof)
% Internal ODE function for timeSimulation.
% State vector: y = [q_reduced; qdot_reduced]

    % extract reduced-order coordinates
    q_red     = y(1:nr);
    q_dot_red = y(nr+1:end);

    % reconstruct physical DOFs for nonlinear force evaluation
    q_phys     = zeros(ndof, 1);
    q_dot_phys = zeros(ndof, 1);
    q_phys(dof)     = T * q_red;
    q_dot_phys(dof) = T * q_dot_red;
   
    % nonlinear bearing forces
    [nonLinBearForce] = nonLinBearingMatrix(Rotor, rotorSpeed, q_phys, q_dot_phys);

    % synchronous rotating excitation: F(t) = Re(F_static * e^{i*Omega*t})
    rot_factor = exp(1i * rotorSpeed * t);
    Force = real(Force_static * rot_factor) + nonLinBearForce + gravityForce;
    ForceR = T_force * Force(dof);

    F = [zeros(nr,1); ForceR];
    ydot = A*y + F;
end

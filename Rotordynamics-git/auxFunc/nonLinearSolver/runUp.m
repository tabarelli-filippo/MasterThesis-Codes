function [time, q, qdot, speed] = runUp(Rotor, t_span, alpha, initPos, initVel, options)
% RUNUP  Computes the time-domain transient response of a rotor during a
%   speed run-up (or coast-down) with a prescribed angular acceleration law.
%
%   The rotor speed varies quadratically in time according to:
%       phi(t)    = alpha(1)*t^2 + alpha(2)*t + alpha(3)    [rad]
%       Omega(t)  = 2*alpha(1)*t + alpha(2)                 [rad/s]
%       dOmega/dt = 2*alpha(1)                              [rad/s^2]
%
%   The nonlinear equations of motion are integrated using ode15s (stiff
%   solver), accounting for:
%     - Speed-dependent gyroscopic and internal damping terms (Omega*C1, Omega*K1)
%     - Nonlinear bearing forces (via nonLinBearingMatrix)
%     - Unbalance and shaft bow excitation (rotated by e^{i*phi(t)})
%     - Angular acceleration contribution to the unbalance force
%
%   Optional modal reduction (Craig-Bampton-style) via reductionMatrix can
%   significantly reduce computation time for large models.
%
% SYNTAX
%   [time, q, qdot, speed] = runUp(Rotor, t_span, alpha, initPos, initVel)
%   [time, q, qdot, speed] = runUp(Rotor, t_span, alpha, initPos, initVel, "nr", nr)
%
% INPUT ARGUMENTS
%   Rotor   - (1x1 struct) Rotor data structure (see rotorMatrix,
%             bearingMatrix, forcing for field definitions)
%   t_span  - (1x2 double) [t_start, t_end] simulation time interval [s]
%   alpha   - (3x1 or 1x3 vector, double) Angular acceleration law
%             coefficients [alpha(1), alpha(2), alpha(3)] such that:
%               phi(t) = alpha(1)*t^2 + alpha(2)*t + alpha(3)
%             For a linear run-up from Omega_0 to Omega_1 in time T:
%               alpha = [(Omega_1-Omega_0)/(2*T), Omega_0, 0]
%   initPos - (ndof x 1 double) Initial nodal displacement vector [m]
%   initVel - (ndof x 1 double) Initial nodal velocity vector [m/s]
%
% NAME-VALUE OPTIONS
%   "nr" - (non-negative integer, default 0) Number of modes to retain in
%          the modal reduction. nr = 0 disables reduction (full model).
%          Reduction is computed at Omega = 0 (static matrices). Choose
%          nr such that the retained bandwidth exceeds the highest speed
%          of interest by a comfortable margin.
%
% OUTPUT ARGUMENTS
%   time  - (ntime x 1 double) Time vector from the ODE solver [s]
%   q     - (ntime x ndof double) Nodal displacement at each time step [m]
%           Rows = time steps, columns = DOFs (zero at constrained DOFs)
%   qdot  - (ntime x ndof double) Nodal velocity at each time step [m/s]
%   speed - (ntime x 1 double) Instantaneous rotor speed Omega(t) [rad/s]
%
% NOTES
%   - ODE solver tolerances: RelTol = 1e-4, AbsTol = 1e-6.
%   - Bearing static matrices are evaluated at Omega = 0 and held fixed;
%     speed-dependent bearing types (e.g., fluid-film type 7) are not
%     updated during the run-up. Use timeSimulation for speed-dependent
%     bearing linearisation at a fixed speed.
%   - ndof = 4 * n_nodes
%
% EXAMPLE
%   ndof = 4 * numel(Rotor.nodes);
%   % Linear run-up from 0 to 1500 rpm in 10 seconds
%   Omega_end = 1500 * pi/30;
%   alpha_law = [Omega_end/20, 0, 0];
%   [t, q, qd, spd] = runUp(Rotor, [0, 10], alpha_law, ...
%                            zeros(ndof,1), zeros(ndof,1), "nr", 8);
%   % Plot speed history
%   figure; plot(t, spd*30/pi); xlabel('t [s]'); ylabel('Speed [rpm]');
%
% SEE ALSO
%   timeSimulation, reductionMatrix, nonLinBearingMatrix, forcing,
%   plotBifurcation, poincareMap, DFT

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

% active (unconstrained) DOFs
rotSpeedref = 0.0;
[~,~,~,zero_dof,~] = bearingMatrix(Rotor,rotSpeedref);
dof = 1:ndof;
dof(zero_dof) = [];
nzero = length(zero_dof);
ncdof = ndof - nzero;

% global matrices
[M0,C0,C1,K0,K1] = rotorMatrix(Rotor);

% bearing matrices at reference speed (held constant during run-up)
rotSpeedref = 0.0; 
[Mb_st, Cb_st, Kb_st, zero_dof, ~] = bearingMatrix(Rotor, rotSpeedref);

M_tot = M0 + Mb_st;
K_tot = K0 + Kb_st;
C_tot = C0 + Cb_st;

%% Optional modal reduction
if nr > 0 && nr < ncdof
   [T, ~] = reductionMatrix(Rotor, rotSpeedref, nr);

   Mr  = T.' * M_tot(dof,dof) * T;
   Kr  = T.' * K_tot(dof,dof) * T;
   Cr  = T.' * C_tot(dof,dof) * T;
   C1r = T.' * C1(dof,dof) * T; 
   K1r = T.' * K1(dof,dof) * T;
   
   q0    = T.' * initPos(dof);
   qdot0 = T.' * initVel(dof);
else
   nr = ncdof;
   T  = eye(ncdof);
   
   Mr  = M_tot(dof,dof);
   Kr  = K_tot(dof,dof);
   Cr  = C_tot(dof,dof);
   C1r = C1(dof,dof);
   K1r = K1(dof,dof);
   
   q0    = initPos(dof);
   qdot0 = initVel(dof);
end

% pre-evaluate constant state-space blocks
A21_st  = -Mr\Kr;
A22_st  = -Mr\Cr;
A21_dyn = -Mr\K1r;   % gyroscopic + internal damping speed-dependent part
A22_dyn = -Mr\C1r;

% forcing
T_force = Mr\T.';
[unForce, bendForce] = forcing(Rotor);
 
initCond = [q0; qdot0];

%% ODE integration
ode_fun = @(t, y) system(t, y, A21_st, A22_st, A21_dyn, A22_dyn, unForce,...
    bendForce, T_force, Rotor, alpha, T, dof, nr, ndof);

ode_opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
[time, y] = ode15s(ode_fun, t_span, initCond, ode_opts);

% instantaneous speed from the prescribed law
speed = 2*alpha(1)*time + alpha(2)*ones(size(time));

% reconstruct full-DOF physical solution
q    = zeros(length(time), ndof);
qdot = zeros(length(time), ndof);

if nr < ncdof
    q(:, dof)    = y(:, 1:nr) * T.';
    qdot(:, dof) = y(:, nr+1:end) * T.';
else
    q(:, dof)    = y(:, 1:nr);
    qdot(:, dof) = y(:, nr+1:end);
end
end

%% -------------------------------------------------------------------------
%  ODE RIGHT-HAND SIDE
% -------------------------------------------------------------------------
function [ydot] = system(t, y, A21_st, A22_st, A21_dyn, A22_dyn, unForce,...
    bendForce, T_force, Rotor, alpha, T, dof, nr, ndof)
% Internal ODE function for runUp.
% State vector: y = [q_reduced; qdot_reduced]

    % angular kinematics
    phi       = alpha(1)*t^2 + alpha(2)*t + alpha(3);
    phidot    = 2*alpha(1)*t + alpha(2);       % Omega(t)
    phidotdot = 2*alpha(1);                    % dOmega/dt

    % speed-dependent state-space blocks
    A21 = A21_st + phidot*A21_dyn;
    A22 = A22_st + phidot*A22_dyn;

    % extract reduced coordinates and velocities
    q_red     = y(1:nr);
    q_dot_red = y(nr+1:end);
    
    % reconstruct physical DOFs for nonlinear force evaluation
    q_phys     = zeros(ndof, 1);
    q_dot_phys = zeros(ndof, 1);
    q_phys(dof)     = T * q_red;
    q_dot_phys(dof) = T * q_dot_red;
    
    % nonlinear bearing forces
    [nonLinBearForce] = nonLinBearingMatrix(Rotor, phidot, q_phys, q_dot_phys);
    
    % rotating unbalance force (acceleration term included for run-up)
    rot_factor = exp(1i * phi);
    ff = bendForce + unForce*(phidot^2 - 1i*phidotdot); 
    Force = real(ff*rot_factor) + nonLinBearForce;
    ForceR = T_force * Force(dof);

    % assemble state-space RHS
    A   = [zeros(nr,nr) eye(nr); A21, A22];
    F   = [zeros(nr,1); ForceR];
    
    ydot = A*y + F;
end

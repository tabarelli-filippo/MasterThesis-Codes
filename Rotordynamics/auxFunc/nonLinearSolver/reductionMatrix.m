function [T, maxPuls] = reductionMatrix(Rotor, rotorSpeed, nr)
% REDUCTIONMATRIX  Computes the modal reduction (Craig-Bampton-style)
%   transformation matrix for a rotor-bearing system at a given speed.
%
%   The reduction is based on the undamped, un-gyroscopic eigenmodes of the
%   system at the specified speed. Gyroscopic and damping contributions are
%   neglected when computing the mode shapes, which makes the approach
%   suitable as a starting-point reduction for weakly gyroscopic systems or
%   as a computationally efficient approximation for time-domain simulation.
%
%   The reduced model has nr degrees of freedom (modal coordinates), down
%   from ncdof physical DOFs (after removing constrained DOFs). The
%   reduction matrix T maps modal coordinates q_r to physical DOFs:
%       q_physical(dof) ≈ T * q_r,   q_r in R^nr
%
% SYNTAX
%   [T, maxPuls] = reductionMatrix(Rotor, rotorSpeed, nr)
%
% INPUT ARGUMENTS
%   Rotor      - (1x1 struct) Rotor data structure (see rotorMatrix,
%                bearingMatrix for field definitions)
%   rotorSpeed - (scalar double) Shaft speed [rad/s] at which the
%                speed-dependent stiffness correction (Omega*K1) and
%                bearing matrices are evaluated. Damping and gyroscopic
%                terms are excluded from the eigenvalue problem.
%   nr         - (scalar double) Number of modes to retain in the reduced
%                model. Must satisfy 1 <= nr < ncdof.
%
% OUTPUT ARGUMENTS
%   T       - (ncdof x nr double) Modal transformation matrix. Columns
%             are the first nr mass-normalised undamped eigenvectors of
%             the constrained system, ordered by ascending natural frequency.
%   maxPuls - (scalar double) Natural frequency of the highest retained
%             mode [Hz]. Provides an estimate of the reduced model's
%             frequency bandwidth: the model is reliable up to ~maxPuls Hz.
%
% NOTES
%   - Constrained DOFs (from bearing types 1 and 2) are removed before
%     computing T; the full-DOF physical vector must be reconstructed by
%     inserting zeros at the constrained positions.
%   - The eigenvalue problem solved is: K * phi = omega^2 * M * phi,
%     where K = K0 + Kb + Omega*K1 and M = M0 + Mb (no C or C1).
%   - For strongly gyroscopic rotors (high speed, large Ip/Id ratio),
%     consider using the full-order model or a gyroscopic-aware reduction.
%   - ndof = 4 * n_nodes;  ncdof = ndof - nzero (constrained DOFs removed)
%
% EXAMPLE
%   [T, f_max] = reductionMatrix(Rotor, 1000, 10);
%   fprintf('Reduced model bandwidth: %.1f Hz\n', f_max);
%   % Use T in timeSimulation via the nr option:
%   [t, q, qdot] = timeSimulation(Rotor, [0,5], 1000, q0, qdot0, "nr", 10);
%
% SEE ALSO
%   timeSimulation, runUp, rotorMatrix, bearingMatrix, charRoots

arguments (Input)
    Rotor (1,1) struct
    rotorSpeed (1,1) double
    nr (1,1) double
end

n_nodes = numel(Rotor.nodes);
ndof = 4 * n_nodes;
dof = 1:ndof;

% global undamped stiffness and mass (gyroscopic and damping excluded)
[M0,~,~,K0,K1] = rotorMatrix(Rotor);
[Mb,~,Kb,zero_dof,~] = bearingMatrix(Rotor,rotorSpeed);

dof(zero_dof) = [];

M = M0 + Mb;
K = K0 + Kb + rotorSpeed*K1;

% solve undamped eigenvalue problem
[Modes, omega_n] = eig(K(dof,dof), M(dof,dof));
[omega_n, isort] = sort(diag(omega_n));
Modes = Modes(:, isort);

% retain the first nr modes
T = Modes(:, 1:nr);
maxPuls = sqrt(omega_n(nr)) / (2*pi);   % highest retained frequency [Hz]

end

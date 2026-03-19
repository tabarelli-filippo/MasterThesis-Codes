function [response] = FRF(Rotor,rotorSpeed)
% FRF  Computes the steady-state forced response of a rotor-bearing system
%   as a function of rotational speed (Frequency Response Function).
%
%   Solves the linear system of equations in the frequency domain at each
%   requested speed, accounting for mass unbalance, slope unbalance, and
%   shaft bow excitation. The equation of motion in the frequency domain is:
%
%       (-Omega^2 * M + i*Omega * C + K) * X = F(Omega)
%
%   where:
%     M = M0 + Mb          (shaft/disk + bearing mass matrices)
%     C = C0 + Cb + Omega*C1 (structural damping + bearing + gyroscopic)
%     K = K0 + Kb + Omega*K1 (elastic + bearing + internal damping)
%     F = bendForce + unForce * Omega^2
%
%   Constrained DOFs (pinned/clamped bearings) are reduced out before
%   solving and restored to zero in the output.
%
% SYNTAX
%   response = FRF(Rotor, rotorSpeed)
%
% INPUT ARGUMENTS
%   Rotor       - (1x1 struct) Rotor data structure containing node, shaft,
%                 disk, bearing, and forcing definitions
%   rotorSpeed  - (vector, double) Array of rotational speeds [rad/s] at
%                 which the response is evaluated
%
% OUTPUT ARGUMENTS
%   response - (ndof x nspeed complex double) Complex displacement response
%              at each DOF and each speed. The amplitude |response| gives
%              the peak displacement; the phase angle(response) gives the
%              phase lag relative to the excitation.
%              DOF ordering: [u1, v1, theta_u1, theta_v1, u2, ...]
%              Units: [m] for translational DOFs, [rad] for rotational DOFs.
%
% NOTES
%   - Gravity loading (forcing type 4) is not included in this function;
%     use a static analysis for gravity-induced deflections.
%   - The response at constrained DOFs is always zero.
%   - ndof = 4 * n_nodes
%   - For orbit plots at a given speed, use:
%       abs(response(:, ispeed))  → amplitude
%       angle(response(:, ispeed)) → phase
%
% EXAMPLE
%   speeds = linspace(10, 1500, 500);   % [rad/s]
%   resp = FRF(Rotor, speeds);
%   % Plot response amplitude at node 3 (x-direction)
%   node = 3; dof_x = 4*node - 3;
%   figure;
%   plot(speeds * 60/(2*pi), abs(resp(dof_x, :)) * 1e6);
%   xlabel('Speed [rpm]'); ylabel('Amplitude [µm]');
%   title('Unbalance response - Node 3, x-direction');
%
% SEE ALSO
%   forcing, rotorMatrix, bearingMatrix, critSpeeds, charRoots

arguments (Input)
    Rotor (1,1) struct
    rotorSpeed (:,:) double {mustBeVector}
end
j = 1i;
n_nodes = numel(Rotor.nodes);
nspeed = length(rotorSpeed);

ndof = 4 * n_nodes;

% define global matrices
[M0,C0,C1,K0,K1] = rotorMatrix(Rotor);
[~,~,~,zero_dof,~] = bearingMatrix(Rotor,rotorSpeed(1));

% sort out zeroed DoF
dof = 1:ndof;
dof(zero_dof) = [];

%% response
[unForce,bendForce] = forcing(Rotor);
response = zeros(ndof,nspeed);
for ispeed = 1:nspeed
   [Mb,Cb,Kb,~] = bearingMatrix(Rotor,rotorSpeed(ispeed));
   M = M0 + Mb;
   K = K0 + Kb + rotorSpeed(ispeed)*K1;
   C = C0 + Cb + rotorSpeed(ispeed)*C1;
   force = bendForce + unForce*(rotorSpeed(ispeed)^2);
   response(dof,ispeed) = (-M(dof,dof)*(rotorSpeed(ispeed)^2)+j*C(dof,dof)*rotorSpeed(ispeed)+K(dof,dof))\force(dof);
end
end

function [response] = FRF(Rotor,rotorSpeed)
%FRF - Frequency Response Function. Evaluates the steady state forced response
%      to unbalance (both mass and slope) and bends
%
%INPUT: Rotor       Structure
%       rotorSpeed  Array of speeds
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

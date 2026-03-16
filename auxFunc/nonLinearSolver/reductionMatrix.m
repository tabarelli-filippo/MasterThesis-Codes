function [T,maxPuls] = reductionMatrix(Rotor,rotorSpeed,nr)
%REDUCTIONMATRIX computes the reduction matrix for a defined rotor speed.
%   Damping and Gyroscopic effects are neglected
%
%INPUT: Rotor       Structure of the rotor
%       rotorSpeed  speed to which reduction is computed
%       nr          number of degrees reduced
%
%OUTPUT:T           Reduction Matrix
%       maxPuls     Highest frequency of the model

arguments (Input)
    Rotor (1,1) struct
    rotorSpeed (1,1) double
    nr (1,1) double
end

n_nodes = numel(Rotor.nodes);
ndof = 4 * n_nodes;
dof = 1:ndof;
%global matrices
[M0,~,~,K0,K1] = rotorMatrix(Rotor);
[Mb,~,Kb,zero_dof,~] = bearingMatrix(Rotor,rotorSpeed);

dof(zero_dof) = [];

M = M0 + Mb;
K = K0 + Kb + rotorSpeed*K1;

[Modes,omega_n] = eig(K(dof,dof),M(dof,dof));
[omega_n,isort] = sort(diag(omega_n));
Modes = Modes(:,isort);
T = Modes(:,1:nr);
maxPuls = sqrt(omega_n(nr))/(2*pi);

end

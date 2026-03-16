function [M0,C0,C1,K0,K1] = rotorMatrix(Rotor)
%ROTORMATRIX computes the global matrices for a given rotor structure
% INPUT: Rotor - Structure
% 
% OUTPUT:M0  mass matrix - global
%        C0  damping matrix - global
%        C1  gyroscopic matrix - global
%        K0  stiffness matrix - global
%        K1  speed dependent contribution to the  stiffness matrix due to the 
%            internal damping - global

n_nodes = Rotor.nodes(end).node;
Shaft = Rotor.shaft;
Disk = Rotor.disk;
n_disks = numel(Disk);

ndof = 4*n_nodes;
% Matrix Initialization
M0 = zeros(ndof, ndof);   % Mass
C0 = zeros(ndof, ndof);   % Damping
C1 = zeros(ndof, ndof);   % Gyroscopic
K0 = zeros(ndof, ndof);   % Stiffness
K1 = zeros(ndof, ndof);   % Speed-dependent Stiffness

%Shaft Matrix Definition
for ii=1:(n_nodes - 1)
    % Element Length
    L = Rotor.nodes(ii+1).coord - Rotor.nodes(ii).coord;
    % Shaft Element Matrix
    if Shaft(ii).type > 0 && Shaft(ii).type < 9.5  
        [M0e,C1e,K0e,K1e,beta] = shaftElement(Shaft(ii),L);
    elseif Shaft(ii).type > 20.5 && Shaft(ii).type < 29.5
        [M0e,C1e,K0e,K1e] = shaftTaperEl(Shaft(ii),L);
    else
        error('Shaft Element at nodes: %d - %d not found\n', Rotor.nodes(ii).node, Rotor.nodes(ii+1).node)
    end
    
    % Dof evaluation
    n1 = Shaft(ii).node1;
    n2 = Shaft(ii).node2;
    dof = [4*n1-3:4*n1 4*n2-3:4*n2];
    % Substitution in the Global Matrix
    M0(dof,dof) = M0(dof,dof) + M0e;
    C0(dof,dof) = C0(dof,dof) + beta*K0e;
    C1(dof,dof) = C1(dof,dof) + C1e;
    K0(dof,dof) = K0(dof,dof) + K0e;
    K1(dof,dof) = K1(dof,dof) +  beta*K1e;
end

% Disk Matrix Definition
for ii = 1:n_disks
    n1 = Disk(ii).node;
    [M0e,C1e] = diskElement(Disk(ii));
    dof = (4*n1-3):4*n1;
    M0(dof,dof) = M0(dof,dof) + M0e;
    C1(dof,dof) = C1(dof,dof) + C1e;
end

end
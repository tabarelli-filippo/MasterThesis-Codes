function [M0,C0,C1,K0,K1] = rotorMatrix(Rotor)
% ROTORMATRIX  Assembles the global finite element matrices of a rotor from
%   shaft element and disk element contributions.
%
%   This is the core assembly function of the FE rotor dynamics solver.
%   It builds the global mass, damping, gyroscopic, and stiffness matrices
%   by looping over all shaft elements and disk elements and adding each
%   element contribution to the appropriate global DOF positions.
%
%   The global equation of motion in the rotating frame is:
%       M0 * x'' + (C0 + Omega*C1) * x' + (K0 + Omega*K1) * x = F
%
%   where:
%     M0        - constant mass matrix (shaft + disks)
%     C0        - structural damping matrix (proportional to shaft K0)
%     C1        - gyroscopic matrix (speed-dependent, proportional to Omega)
%     K0        - elastic stiffness matrix
%     K1        - speed-dependent stiffness correction due to internal
%                 (material) damping (skew-symmetric destabilizing term)
%
%   Bearing contributions (Mb, Cb, Kb) are assembled separately in
%   bearingMatrix and added to these matrices before solving.
%
% SYNTAX
%   [M0, C0, C1, K0, K1] = rotorMatrix(Rotor)
%
% INPUT ARGUMENTS
%   Rotor - (struct) Rotor data structure. Required fields:
%             .nodes  - array of node structs, each with:
%                         .node   (integer node index)
%                         .coord  (axial coordinate z [m])
%             .shaft  - array of shaft element structs, each with:
%                         .type   (integer: 1-9 for uniform, 21-29 tapered)
%                         .node1, .node2 (endpoint node indices)
%                         plus geometry and material fields (see
%                         shaftElement, shaftTaperEl)
%             .disk   - array of disk structs, each with:
%                         .node   (attachment node index)
%                         plus type-specific fields (see diskElement)
%
% OUTPUT ARGUMENTS
%   M0 - (ndof x ndof double) Global mass matrix [kg, kg·m²]
%   C0 - (ndof x ndof double) Global structural damping matrix [N·s/m].
%        Built as C0 = sum_e (beta_e * K0e), where beta_e is the element
%        material damping coefficient [s]
%   C1 - (ndof x ndof double) Global gyroscopic matrix [kg·m²].
%        The speed-dependent term is Omega * C1 in the EOM.
%   K0 - (ndof x ndof double) Global elastic stiffness matrix [N/m]
%   K1 - (ndof x ndof double) Global speed-dependent stiffness matrix.
%        Represents the destabilizing effect of internal damping:
%        K1 = sum_e (beta_e * K1e). The speed-dependent term is Omega*K1.
%
% NOTES
%   - Shaft element types 1-9 are handled by shaftElement (uniform cross-
%     section, Euler-Bernoulli or Timoshenko beam).
%   - Shaft element types 21-29 are handled by shaftTaperEl (linearly
%     tapered cross-section).
%   - ndof = 4 * n_nodes, with DOF ordering [u, v, theta_u, theta_v]
%     at each node.
%   - Bearing matrices must be added separately via bearingMatrix before
%     solving the eigenvalue problem or frequency response.
%
% EXAMPLE
%   [M0, C0, C1, K0, K1] = rotorMatrix(Rotor);
%   % Evaluate system matrices at Omega = 500 rad/s (no bearings)
%   Omega = 500;
%   M_sys = M0;
%   C_sys = C0 + Omega * C1;
%   K_sys = K0 + Omega * K1;
%
% SEE ALSO
%   shaftElement, shaftTaperEl, diskElement, bearingMatrix, charRoots, FRF

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

function [nu_static,nu_fatigue,N_life] = stressAnalysis(Rotor,response,sigma_y)
% STRESSANALYSIS  Performs a static stress analysis on the rotor shaft at
%   each element, evaluating safety factors against yielding under combined
%   bending, axial, and torsional loading.
%
%   For each shaft element, the following stress components are computed
%   at both end nodes (nodes i and i+1):
%     - Bending stress sigma_F: due to lateral curvature of the shaft,
%       evaluated from the deflected shape via the element curvature
%       (second derivative of the Hermite interpolation).
%     - Axial stress sigma_N: due to a constant axial force along the
%       element (tension or compression).
%     - Torsional shear stress tau_T: due to a constant torque along
%       the element.
%
%   The Von Mises equivalent stress is then used to compute the static
%   safety factor at each node:
%       sigma_VM = sqrt(sigma^2 + 3*tau^2)
%       nu_static = sigma_y / sigma_VM
%
% SYNTAX
%   [nu_static, nu_fatigue, N_life] = stressAnalysis(Rotor, response, sigma_y)
%
% INPUT ARGUMENTS
%   Rotor    - (1x1 struct) Rotor data structure containing:
%                .shaft  - array of shaft element structs with geometry
%                          and material fields (.E, .G, .d_ext1, .d_int1,
%                          .d_ext2, .d_int2, .node1); optionally .AxialForce
%                          and .Torque
%                .nodes  - node coordinate array
%   response - (ndof x 1 double or complex double) Nodal displacement
%              vector from a static or forced response analysis (e.g. from
%              FRF at a given speed). Only the real part is used.
%              DOF ordering: [u1, v1, theta_u1, theta_v1, u2, ...]
%   sigma_y  - (positive double, scalar) Material yield stress [Pa]
%
% OUTPUT ARGUMENTS
%   nu_static  - ((n_shaft+1) x 1 double) Static safety factor at each
%                node (shaft nodes only). nu_static > 1 indicates no
%                yielding; nu_static < 1 indicates yielding.
%   nu_fatigue - Reserved for future implementation (currently empty [])
%   N_life     - Reserved for future implementation (currently empty [])
%
% NOTES
%   - Bending stress is computed from the curvature of the deflected shaft
%     using the Euler-Bernoulli curvature-moment relationship:
%       sigma_F = -E * (d/2) * d²w/dz²
%     evaluated at both ends of each element using the Hermite shape
%     function second derivatives.
%   - At each node, the maximum stress among the contributions from the
%     adjacent elements is used.
%   - Disk stress analysis is not yet implemented.
%   - ndof = 4 * n_nodes
%
% EXAMPLE
%   speeds = linspace(100, 2000, 200);
%   resp = FRF(Rotor, speeds);
%   % Stress analysis at first critical speed response
%   [~, idx] = max(abs(resp(dof_x, :)));
%   [nu_s, ~, ~] = stressAnalysis(Rotor, real(resp(:, idx)), 250e6);
%   fprintf('Minimum static safety factor: %.2f\n', min(nu_s));
%
% SEE ALSO
%   FRF, charRoots, rotorMatrix, shaftElement

arguments (Input)
    Rotor (1,1) struct
    response (:,:) double 
    sigma_y (1,1) double {mustBePositive}
end
Shaft = Rotor.shaft;
%% Shaft Stress Analysis
n_shaft = numel(Shaft);
sigma_F = zeros(n_shaft,4); sigma_N = zeros(n_shaft, 2); tau_T = zeros(n_shaft, 2);

if isfield(Shaft, 'AxialForce')
    AxialForce = [Shaft.AxialForce];
else
    AxialForce = zeros(n_shaft,1);
end

if isfield(Shaft, 'Torque')
    Torque = [Shaft.Torque];
else
    Torque = zeros(n_shaft,1);
end

for ii = 1:n_shaft
    % bending stress from lateral curvature (Euler-Bernoulli)
    Le = Rotor.nodes(ii +1).coord - Rotor.nodes(ii).coord;
    E = Shaft(ii).E;
    d_ext1 = Shaft(ii).d_ext1;

    if isempty(Shaft(ii).d_ext2)
        d_ext2 = d_ext1;
    else
        d_ext2 = Shaft(ii).d_ext2;
    end  
    inode = Shaft(ii).node1;

    uz1_curv = [-6, -4*Le, 6, -2*Le] * real(response([4*inode-3 4*inode 4*inode+1 4*inode+4])) / Le^2;
    uz2_curv = [6, 2*Le, -6, 4*Le] * real(response([4*inode-3 4*inode 4*inode+1 4*inode+4])) / Le^2;
    vz1_curv = [-6, 4*Le, 6, 2*Le] * real(response([4*inode-2 4*inode-1 4*inode+2 4*inode+3],1)) / Le^2;
    vz2_curv = [6, -2*Le, -6, -4*Le] * real(response([4*inode-2 4*inode-1 4*inode+2 4*inode+3],1)) / Le^2;  
    
    sigma_F(ii,1) = -0.5*E*d_ext1*sqrt(uz1_curv^2 + vz1_curv^2);
    sigma_F(ii,2) = -0.5*E*d_ext2*sqrt(uz2_curv^2 + vz2_curv^2);

    % axial normal stress
    d_int1 = Shaft(ii).d_int1;
    if isempty(Shaft(ii).d_int2)
        d_int2 = d_int1;
    else
        d_int2 = Shaft(ii).d_int2;
    end
    d_ext = [d_ext1,d_ext2]; d_int = [d_int1,d_int2];
    A = pi/4 * (d_ext.^2 - d_int.^2);
    sigma_N(ii,:) = AxialForce(ii)./A;
    
    % torsional shear stress
    G = Shaft(ii).G;
    J = pi/32 * (d_ext.^4 - d_int.^4);
    tau_T(ii,:) = [Torque(ii) Torque(ii)]/J.*d_ext/2;
end
%reshaping
sigma_F_node = [[sigma_F(:,1);0],[0;sigma_F(:,2)]]; % nshaft x 2
sigma_N_node = [[sigma_N(:,1);0],[0;sigma_N(:,2)]]; % nshaft x 2
tau_T_node = [[tau_T(:,1);0],[0;tau_T(:,2)]];

% max lateral stress for each node
idx_row = (1:n_shaft+1)';

sigma_node = sigma_F_node + sigma_N_node; % nshaft x 2

[~,idx_col] = max(abs(sigma_node),[],2);
idx_lin= sub2ind(size(sigma_node), idx_row, idx_col);
sigma = sigma_node(idx_lin); % nshaft x 1
 
[~,idx_col] = max(abs(tau_T_node),[],2);
idx_lin= sub2ind(size(tau_T_node), idx_row, idx_col);
tau_T_max = tau_T_node(idx_lin); % nshaft x 1

tau = tau_T_max;

% Von Mises Criterion
VM_stress = sqrt(sigma.^2 + 3*tau.^2);
nu_static = (sigma_y* ones(n_shaft+1,1)) ./ VM_stress;

%% Fatigue life analysis - reserved for future implementation
nu_fatigue =[];
N_life =[];

%% Disk Stress Analysis - reserved for future implementation

end

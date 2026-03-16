function [nu_static,nu_fatigue,N_life] = stressAnalysis(Rotor,response,sigma_y)
%STRESSANALISYS computes stress analysis at each node to check material
%failure
%
%INPUT: Rotor       Structure of the rotor
%       response    Steady response - deformated structure
%       sigma_y     Yielding stress [Pa]

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
    % sigma lateral stress
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

    % axial sigma stress
    d_int1 = Shaft(ii).d_int1;
    if isempty(Shaft(ii).d_int2)
        d_int2 = d_int1;
    else
        d_int2 = Shaft(ii).d_int2;
    end
    d_ext = [d_ext1,d_ext2]; d_int = [d_int1,d_int2];
    A = pi/4 * (d_ext.^2 - d_int.^2);
    sigma_N(ii,:) = AxialForce(ii)./A;
    
    %torque tau stress
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

%% implementare la parte di vita a fatica
nu_fatigue =[];
N_life =[];

%% Disk Stress Analysis

end

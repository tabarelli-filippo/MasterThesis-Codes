function [M_red, K_red] = ROM(M, K, Psi_disk, B, N)
%ASSEMBLE_FULL_ROTOR full matrix from sector
%   INPUT:
%   - M, K:      Full global matrices [N*Dofs x N*Dofs]
%   - Psi_disk:  Global Disk mode Matrix [N*Dofs x N_disk_modes]
%   - B:         Blade mode [Dofs_pala x N_blade_modes] 
%   - Kd_vals:   Disk eigenvalues
%   - N:         number of Blades 
%
%   OUTPUT:
%   - M_red, K_red: Matrici ridotte [(N_disk + N) x (N_disk + N)]

Psi_blade = sparse(kron(eye(N), B));

T = [Psi_disk, Psi_blade];

% M_red: T' * M * T
M_red = T' * M * T;
% K_red: T' * M * T
K_red = T' * K * T;

end
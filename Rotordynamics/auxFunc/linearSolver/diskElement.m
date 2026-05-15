function [M0e,C1e] = diskElement(Disk)
% DISKELEMENT  Computes the element mass and gyroscopic matrices for a
%   rigid disk attached to the rotor.
%
%   Two definition modes are supported:
%     Type 1 - Geometric definition: disk inertia properties are derived
%              from external diameter, internal diameter, thickness, and
%              material density.
%     Type 2 - Direct definition: mass and inertia values are provided
%              explicitly.
%
%   The disk is modelled as a rigid body with 4 DOFs at its node:
%   [u, v, theta_u, theta_v], where u and v are lateral displacements
%   and theta_u, theta_v are the corresponding bending rotations.
%
% SYNTAX
%   [M0e, C1e] = diskElement(Disk)
%
% INPUT ARGUMENTS
%   Disk - (struct) Disk definition structure. Required fields depend on
%          the disk type:
%
%          Type 1 (geometric):
%            .type   = 1
%            .D_ext  - External diameter [m]
%            .D_int  - Internal diameter [m]
%            .thick  - Axial thickness [m]
%            .rho    - Material density [kg/m³]
%
%          Type 2 (direct inertia input):
%            .type   = 2
%            .mass   - Total disk mass [kg]
%            .Ip     - Polar moment of inertia [kg·m²]
%            .Id     - Diametral moment of inertia [kg·m²]
%
% OUTPUT ARGUMENTS
%   M0e - (4x4 double) Element mass matrix in local DOF ordering
%         [u, v, theta_u, theta_v]. Diagonal entries are m, m, Id, Id.
%   C1e - (4x4 double) Element gyroscopic matrix. The gyroscopic coupling
%         terms involve the polar moment of inertia Ip; the matrix is
%         skew-symmetric and proportional to the shaft speed when assembled
%         into the global system via rotorMatrix.
%
% THEORY
%   For a thin rigid disk:
%     m  = rho * pi/4 * thick * (D_ext^2 - D_int^2)
%     Ip = rho * pi/32 * thick * (D_ext^4 - D_int^4)
%     Id = Ip/2 + m * thick^2 / 12
%
%   The gyroscopic matrix contribution to the equations of motion is
%   Omega * C1e, where Omega is the shaft rotational speed.
%
% EXAMPLE
%   Disk.type  = 1;
%   Disk.D_ext = 0.3;    % [m]
%   Disk.D_int = 0.05;   % [m]
%   Disk.thick = 0.04;   % [m]
%   Disk.rho   = 7800;   % [kg/m^3]
%   [M0e, C1e] = diskElement(Disk);
%
% SEE ALSO
%   rotorMatrix, shaftElement, shaftTaperEl

type = Disk.type;
switch type
    case 1
        D_ext = Disk.D_ext;
        D_int = Disk.D_int;
        thick = Disk.thick;
        rho = Disk.rho;
        
        m = rho * pi / 4 * thick * (D_ext^2 - D_int^2);
        Ip = rho * pi / 32 * thick * (D_ext^4 - D_int^4);
        Id = Ip / 2 + m * (thick^2) / 12;
    case 2
        m = Disk.mass;
        Ip = Disk.Ip;
        Id = Disk.Id;
end
    M0e = diag([m m Id Id]);
    C1e = zeros(4,4);
    C1e(3,4) = Ip;
    C1e(4,3) = -Ip;
end

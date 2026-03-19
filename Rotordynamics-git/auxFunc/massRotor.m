function mass = massRotor(Rotor)
% MASSROTOR  Computes the total mass of the rotor from its geometric and
%   material description, summing contributions from shaft elements and
%   disks.
%
%   Shaft element masses are computed analytically for linearly tapered
%   hollow cylinders using the exact formula for a truncated cone:
%       m_e = rho * pi/12 * L * (d_ext2^2 + d_ext1*d_ext2 + d_ext1^2
%                               - d_int2^2 - d_int1*d_int2 - d_int1^2)
%   For uniform elements (d1 = d2), this reduces to the standard cylinder
%   formula m = rho * A * L.
%
%   Disk masses are either computed from geometry (type 1) or read directly
%   from the struct (type 2).
%
% SYNTAX
%   mass = massRotor(Rotor)
%
% INPUT ARGUMENTS
%   Rotor - (1x1 struct) Rotor data structure. Required fields:
%             .nodes  - node array with .coord (axial coordinates [m])
%             .shaft  - array of shaft element structs with:
%                         .d_ext1, .d_ext2  external diameters [m]
%                         .d_int1, .d_int2  internal diameters [m]
%                         .rho              density [kg/m³]
%             .disk   - array of disk structs:
%                         Type 1: .D_ext, .D_int [m], .thick [m], .rho [kg/m³]
%                         Type 2: .mass [kg]
%
% OUTPUT ARGUMENTS
%   mass - (scalar double) Total rotor mass [kg] = shaft mass + disk mass
%
% EXAMPLE
%   m_total = massRotor(Rotor);
%   fprintf('Total rotor mass: %.3f kg\n', m_total);
%
% SEE ALSO
%   diskElement, shaftElement, meshGenerator, figureRotor

arguments (Input)
    Rotor (1,1) struct
end

Nodes = Rotor.nodes; Shaft = Rotor.shaft; Disk = Rotor.disk;

nodes_coord = [Nodes.coord];

%% Shaft mass (truncated cone formula per element)
d_int1 = [Shaft.d_int1];    d_ext1 = [Shaft.d_ext1];
d_int2 = [Shaft.d_int2];    d_ext2 = [Shaft.d_ext2];
rhoShaft = [Shaft.rho];

lengthEl = nodes_coord(2:end) - nodes_coord(1:end-1);
massShaftEl = pi/12 * rhoShaft .* lengthEl .* ...
    (d_ext2.^2 + d_ext1.*d_ext2 + d_ext1.^2 ...   % external cone
   - d_int2.^2 - d_int1.*d_int2 - d_int1.^2);     % internal cone (hollow)

%% Disk mass
disk_nodes = [Disk.node];
numDisks = length(disk_nodes); 
type = [Disk.type];
massDisk = zeros(numDisks,1);

for ii = 1:numDisks
    switch type(ii)
        case 1  % geometric definition
            D_ext  = Disk(ii).D_ext;
            D_int  = Disk(ii).D_int;
            thick  = Disk(ii).thick;
            rhoDisk = Disk(ii).rho;
            massDisk(ii) = rhoDisk * pi/4 * (D_ext^2 - D_int^2) * thick;
        case 2  % direct mass input
            massDisk(ii) = Disk(ii).mass;
    end
end

%% Total mass
mass = sum(massShaftEl) + sum(massDisk);
end

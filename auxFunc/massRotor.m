function mass = massRotor(Rotor)
%MASSTROTOR evaluates Rotor's mass base on description
arguments (Input)
    Rotor (1,1) struct
end
%%
Nodes = Rotor.nodes; Shaft = Rotor.shaft; Disk = Rotor.disk;

nodes_coord = [Nodes.coord];
%% Shaft
d_int1 = [Shaft.d_int1];    d_ext1 = [Shaft.d_ext1];
d_int2 = [Shaft.d_int2];    d_ext2 = [Shaft.d_ext2];
rhoShaft = [Shaft.rho];

% element mass
lengthEl = nodes_coord(2:end) - nodes_coord(1:end-1);
massShaftEl = pi/12*rhoShaft.*lengthEl.*(d_ext2.^2 + d_ext1.*d_ext2 + d_ext2.^2 ... %external volume - internal volume
    - d_int2.^2 - d_int1.*d_int2 - d_int1.^2);

%% Disks
disk_nodes = [Disk.node];
numDisks = length(disk_nodes); 
type = [Disk.type];
massDisk = zeros(numDisks,1);

for ii = 1:numDisks
    switch type(ii)
        case 1
            D_ext = Disk(ii).D_ext;
            D_int = Disk(ii).D_int;
            thick = Disk(ii).thick;
            rhoDisk = Disk(ii).rho;
            massDisk(ii) = rhoDisk*pi/4*(D_ext^2 - D_int^2)*thick;
        case 2
            massDisk(ii) = Disk(ii).mass;
    end
end

%% Total mass
mass = sum(massShaftEl) + sum(massDisk);
end
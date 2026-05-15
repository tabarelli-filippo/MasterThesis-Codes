% MAIN  Modal analysis of a turbocharger rotor modelled from geometry with
%   the meshGenerator utility and two pinned bearings.
%
%   Builds a refined FE mesh for a small turbocharger shaft from measured
%   shoulder coordinates and disk/bearing positions. Performs Campbell
%   diagram analysis and root locus for the first seven modes.
%
%   ROTOR DESCRIPTION
%     Shaft  : 13 coarse sections, each refined with 2 intermediate nodes
%              (Timoshenko type 3, except element 13 which is tapered type 23)
%              Total length ≈ 234 mm; diameters range from 12 to 68.6 mm
%     Disks  : 2 rigid disks (type 2 — direct mass/inertia input):
%                Compressor: m = 2.3 kg, Id = 3.5e-6 kg·m², Ip = 4.0e-6 kg·m²
%                            at z = 50 mm
%                Turbine:    m = 3.65 kg, Id = 4.0e-6 kg·m², Ip = 8.0e-6 kg·m²
%                            at z ≈ 215 mm
%     Bearings: 2 pinned (type 1) at effective positions computed from
%               ball bearing contact geometry (mean diameter d_m, contact
%               angle alpha = 20°), nominally at z ≈ 101 mm and z ≈ 169 mm
%
%   ANALYSIS PERFORMED
%     1. Mesh generation with meshGenerator (13 coarse elements,
%        2 intermediate nodes each, element 13 tapered)
%     2. Disk and bearing assembly
%     3. Rotor schematic and total mass
%     4. Campbell diagram for 7 sorted modes, 1X and 2X lines
%        (charRoots, sortModesMAC, plotCampbell)
%     5. Mode shape visualisation at the 5th speed step (plotMode)
%     6. Root locus for the first mode (plotRootLocus)
%
% SEE ALSO
%   meshGenerator, charRoots, sortModesMAC, plotCampbell, plotMode,
%   plotRootLocus, massRotor, figureRotor

clear
close all
clc
addpath(genpath('auxFunc/'))

%% Shaft shoulder coordinates and diameter table
coord_should = [0 20 34 42 52.9 93.9 108.3 160.3 176.5 184.7 185.7 187.7 193.9 233.3]'*1e-3;
d_Ext1Shaft  = [6 7 6.2 9.5 9.825 10 9.825 10 12.7 19.3 34.3 34.3 18.5]' * 2e-3;
d_Ext2Shaft  = [6 7 6.2 9.5 9.825 10 9.825 10 12.7 21.3 34.3 34.3 18.5]' * 2e-3;
d_IntShaft   = zeros(size(d_Ext1Shaft));

rho     = 7800;
E       = 200e9;
poisson = 0.285;
G       = E / 2 / (1 + poisson);

%% Disk parameters
z_compr      = 50e-3;
m_compr      = 2300e-3; Id_compr = 3.5e6*1e-9; Ip_compr = 4e6*1e-9;
d_ext_shaft_comp = 9.5*2e-3; d_int_shaft_comp = 0;

z_turb       = 214.71e-3;
m_turb       = 3650e-3; Id_turb = 4e6*1e-9; Ip_turb = 8e6*1e-9;
d_extShaft_turb  = 18.5*2e-3; d_intShaft_turb = 0;

z_disks = [z_compr; z_turb];

%% Bearing positions (from ball bearing contact geometry)
d_m_bearing = 0.3e-3;
alpha       = 20*pi/180;
delta       = 0.5*d_m_bearing*tan(alpha);
z_bearings  = [(93.9+108.3)/2-delta; (160.3+176.5)/2+delta]*1e-3;
d_extShaft_bear = 10*2e-3; d_intShaft_bear = 0;

%% Reorder coordinates and build diameter table
[z_coords, idx_sort] = sort([coord_should; z_disks; z_bearings]);
numEl    = length(z_coords) - 1;
zero_vec = zeros(numEl, 1);

d_extShaft_augm1 = [d_Ext1Shaft; d_ext_shaft_comp; d_extShaft_turb; d_extShaft_bear; d_extShaft_bear];
d_extShaft_augm2 = [d_Ext2Shaft; d_ext_shaft_comp; d_extShaft_turb; d_extShaft_bear; d_extShaft_bear];
d_intShaft_augm  = [d_IntShaft;  d_int_shaft_comp;  d_intShaft_turb;  d_intShaft_bear;  d_intShaft_bear];

d_shaftEl        = [zero_vec, d_extShaft_augm1, zero_vec, d_extShaft_augm2];
d_shaftEl        = d_shaftEl(idx_sort(2:end)-1, :);

%% Mesh parameters
n_intNodes      = 2*ones(numEl,1);   % 2 intermediate nodes per element
isGraded        = zero_vec;          % uniform spacing
shaftElType     = 3*ones(numEl,1);   % Timoshenko (Hutchinson)
shaftElType(13) = 23;                % element 13: tapered Timoshenko
shaftElProperties = [rho*ones(numEl,1), E*ones(numEl,1), poisson*ones(numEl,1)];

%% Generate mesh
rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, shaftElType, shaftElProperties);

%% Disk and bearing assembly
z_nodes_mesh = [rotor.nodes.node];
z_coord_mesh = [rotor.nodes.coord];

idx_compr = find(z_coord_mesh == z_compr, 1);
idx_turb  = find(z_coord_mesh == z_turb, 1);
rotor.disk(1) = struct('type',2,'node',z_nodes_mesh(idx_compr),'mass',m_compr,'Id',Id_compr,'Ip',Ip_compr);
rotor.disk(2) = struct('type',2,'node',z_nodes_mesh(idx_turb), 'mass',m_turb, 'Id',Id_turb, 'Ip',Ip_turb);

idx_bear1 = find(z_coord_mesh == z_bearings(1), 1);
idx_bear2 = find(z_coord_mesh == z_bearings(2), 1);
rotor.bearing(1) = struct('type',1,'node',z_nodes_mesh(idx_bear1));
rotor.bearing(2) = struct('type',1,'node',z_nodes_mesh(idx_bear2));

%% Rotor schematic and mass
figureRotor(rotor);
mass = massRotor(rotor);
fprintf('Total rotor mass: %.3f kg\n', mass);

%% Campbell diagram
rotorSpeed = (0:1000:10000)*pi/30;

[natPuls, Mode, kappa] = charRoots(rotor, rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

NatPuls_sort0 = NatPuls_sort(1:7,:);
kappa_sort0   = kappa_sort(:,1:7,:);

NX = [1, 2]; isDamped = true;
plotCampbell(rotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort0);

%% Mode shapes at the 5th speed step
figure('Name','Modes','NumberTitle','off')
subplot(1,3,1); plotMode(rotor,Mode_sort(:,1,5),'number',1,'frequency',abs(imag(NatPuls_sort0(1,5)))/2/pi,'rotorSpeed',rotorSpeed(5)*30/pi);
subplot(1,3,2); plotMode(rotor,Mode_sort(:,3,5),'number',3,'frequency',abs(imag(NatPuls_sort0(3,5)))/2/pi,'rotorSpeed',rotorSpeed(5)*30/pi);
subplot(1,3,3); plotMode(rotor,Mode_sort(:,5,5),'number',5,'frequency',abs(imag(NatPuls_sort0(5,5)))/2/pi,'rotorSpeed',rotorSpeed(5)*30/pi);

%% Root locus (first mode)
plotRootLocus(rotorSpeed*30/pi, NatPuls_sort0(1,:));

% VALID_1  Model validation against published rotor dynamics data.
%   Reference: doi:10.1016/j.measurement.2012.01.032
%
%   Validates the FE rotor model against the turbocharger benchmark from the
%   reference paper. The paper uses a lumped (not consistent) mass matrix for
%   the shaft, which produces lower resonance frequencies. This script uses
%   the correct consistent mass formulation, resulting in higher natural
%   frequency values. The paper reports only forward-whirl modes (backward
%   modes plotted as negative); all critical speeds are computed here without
%   FW/BW distinction.
%
%   ROTOR DESCRIPTION
%     Shaft : 10 elements, total length ≈ 453.5 mm
%             hollow sections (elements 1 and 4–6) and solid sections
%             Mixed densities: rho = 7810 kg/m³ for shaft material,
%             rho = 0 for structural elements representing non-metallic parts
%     Disks : 2 rigid disks (type 2, direct inertia input):
%               Compressor: m = 20.81 kg at z = 98.3 mm
%               Turbine:    m = 18.2 kg at z = 363.5 mm
%               Compressor unbalance: epsilon = 1e-6 m
%     Bearings: 2 pinned (type 1) at z = 0 and z = 453.5 mm
%     Elements: Euler-Bernoulli (type 1), no intermediate refinement
%
%   ANALYSIS PERFORMED
%     1. Mesh generation (meshGenerator, no refinement)
%     2. Disk and bearing assembly
%     3. Rotor schematic (figureRotor)
%     4. Campbell diagram for 12 modes, 1X line
%        (charRoots, sortModesMAC, plotCampbell)
%     5. Critical speeds (iterative method 2, 9 critical speeds)
%        and mode shape visualisation for selected modes (plotMode)
%     6. Unbalance frequency response at nodes 3.1 and 3.2 (FRF, plotFRF)
%     7. Deflected shape at 2nd, 4th, and 6th critical speeds
%        (plotDisplacement)
%
% SEE ALSO
%   valid_2, valid_3, valid_4, meshGenerator, charRoots, critSpeeds, FRF

clear
close all
addpath(genpath('auxFunc/'))

%% Shaft geometry
z_coords   = [0 30 98.3 132.5 162.5 250.5 338.5 363.5 393.5 433.5 453.5]'*1e-3;
d_ExtShaft = [30 120 120 63 68 68 150 150 40 30]'*1e-3;
d_IntShaft = [15 0 0 55 60 60 0 0 30 20]'*1e-3;
numEl = length(d_IntShaft);

rho = [7810 0 0 7810 7810 7810 0 0 7810 7810]';
E   = 210e9;

%% Disk parameters
z_compr  = 98.3e-3;  m_compr = 20.81; Id_compr = 0.174; Ip_compr = 0.285;
z_turb   = 363.5e-3; m_turb  = 18.2;  Id_turb  = 0.142; Ip_turb  = 0.269;
z_disks  = [z_compr; z_turb];
epsilon  = 1e-6;     % compressor unbalance eccentricity [m]

%% Bearing positions
z_bearings = [0; 453.5e-3];

%% Mesh (no refinement, Euler-Bernoulli)
zero_vec          = zeros(numEl,1);
n_intNodes        = zero_vec;
isGraded          = zero_vec;
shaftElType       = ones(numEl,1);
shaftElProperties = [rho, E*ones(numEl,1), zero_vec];
d_shaftEl         = [d_IntShaft, d_ExtShaft, d_IntShaft, d_ExtShaft];

rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, shaftElType, shaftElProperties);

%% Disk and bearing assembly
z_nodes_mesh = [rotor.nodes.node];
z_coord_mesh = [rotor.nodes.coord];

idx_compr = find(z_coord_mesh == z_compr, 1);
idx_turb  = find(z_coord_mesh == z_turb, 1);
rotor.disk(1) = struct('type',2,'node',z_nodes_mesh(idx_compr),'mass',m_compr,'Id',Id_compr,'Ip',Ip_compr);
rotor.disk(2) = struct('type',2,'node',z_nodes_mesh(idx_turb), 'mass',m_turb, 'Id',Id_turb, 'Ip',Ip_turb);
rotor.disk(1).epsilon = epsilon;
rotor.disk(2).epsilon = 0;
rotor.forcing(1).type = 1;

idx_bear1 = find(z_coord_mesh == z_bearings(1), 1);
idx_bear2 = find(z_coord_mesh == z_bearings(2), 1);
rotor.bearing(1) = struct('type',1,'node',z_nodes_mesh(idx_bear1));
rotor.bearing(2) = struct('type',1,'node',z_nodes_mesh(idx_bear2));

%% Rotor schematic
figureRotor(rotor);

%% Campbell diagram
rotorSpeed = (0:1e3:1e5);
[natPuls, Mode, kappa] = charRoots(rotor, rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

selectModes   = 1:12;
NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0   = kappa_sort(:,selectModes,:);
NX = 1; isDamped = false;
plotCampbell(rotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort0);

%% Critical speeds and mode shapes
method = 2;
[criticalSpeeds, modeCrits] = critSpeeds(rotor, NX, isDamped, method, 'num_crit', 9);

figure('Name','Modes','NumberTitle','off')
selectModes = [2, 4, 8];
for jj = 1:length(selectModes)
    subplot(1,3,jj)
    kk = selectModes(jj);
    plotMode(rotor, modeCrits(:,kk), 'number', kk, ...
        'frequency', abs(imag(NatPuls_sort0(1,kk)))/2/pi, ...
        'rotorSpeed', criticalSpeeds(kk)*30/pi);
    fprintf('Critical Speed %d: %.3f rpm\n', kk, criticalSpeeds(kk)*30/pi);
end

%% Unbalance response
rotorSpeed = 0:10:6000;
response   = FRF(rotor, rotorSpeed);
plotFRF(rotorSpeed*30/pi, response, [3.1 3.2])

%% Deflected shape at selected critical speeds
idx_crit = find(rotorSpeed == 2310, 1);
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

idx_crit = find(rotorSpeed == 3750, 1);
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

idx_crit = find(rotorSpeed == 5640, 1);
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

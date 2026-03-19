% VALID_4  Model validation against published aero-engine rotor data.
%   Reference: 10.2514/1.B34104
%
%   Validates the FE model against a multi-stage aero-engine rotor benchmark
%   with six bladed disks, a hollow shaft, and two speed-dependent bearings
%   (type 13). Includes a complete run-up simulation from 0 to 500 Hz using
%   the runUp integrator with modal reduction.
%
%   ROTOR DESCRIPTION
%     Shaft  : 15 elements, total length = 780 mm
%              hollow central section (elements 4–10: d_int = 58 mm)
%              outer sections solid; d_ext ranges from 15 mm to 80 mm
%              Material: steel (rho = 7850 kg/m³, E = 205 GPa)
%              Element 10 refined with 6 intermediate nodes (cosine grading off)
%     Disks  : 6 bladed disks (type 1, geometric) at z = elements 5–9 and 15:
%              D_ext = 190–250 mm; mixed aluminum, titanium, steel densities
%              Unbalance: disk 6 only (m0 = 1e-4 kg / 0.125 m)
%     Bearings: 2 speed-dependent (type 13):
%               Front: kxx = 13.7e6 N/m (constant), cxx varies 5.1–39.3 N·s/m
%               Rear:  kxx varies 2.94–12.5e6 N/m, cxx varies 2441–9720 N·s/m
%               Both: diagonal, no cross-coupling
%
%   ANALYSIS PERFORMED
%     1. Mesh generation (Timoshenko type 3, element 10 refined 6 nodes)
%     2. Disk and bearing assembly; rotor schematic
%     3. Campbell diagram for 6 sorted modes, 1X and 0.5X lines (isDamped)
%     4. Critical speeds (iterative method 2, 6 speeds) with mode shapes
%     5. Stability: root locus for FW modes 1 and 2 (plotRootLocus)
%     6. Logarithmic decrement (LogDecrement, plotLogDecrement)
%     7. Unbalance response at nodes 5 and 21 (FRF, plotFRF)
%     8. Deflected shape at 4th critical speed (plotDisplacement)
%     9. Run-up simulation 0→500 Hz in 10 s with 12-mode reduction
%        (runUp, nr = 12); time history and response vs speed plotted
%        at node 21
%
% SEE ALSO
%   valid_1, valid_2, valid_3, charRoots, critSpeeds, runUp, plotLogDecrement

clear
close all
addpath(genpath('auxFunc/'))

%% Shaft geometry
z_coords = [0;0.04;0.075;0.087;0.1;0.157;0.21;0.253;0.292;0.36;...
            0.63;0.642;0.7;0.735;0.765;0.78];
numEl = length(z_coords)-1;

d_ExtShaft = [0.015;0.015;0.040;0.040;0.040;0.040;0.040;...
              0.040;0.040;0.034;0.034;0.015;0.015;0.0275;0.0275]*2;
d_IntShaft = [0;0;0;0.029;0.029;0.029;0.029;0.029;0.029;...
              0.029;0;0;0;0;0]*2;

rho     = 7850*ones(numEl,1);
E       = 2.05e11*ones(numEl,1);
poisson = 0.29*ones(numEl,1);
beta    = zeros(numEl,1);

%% Disk parameters
z_disks        = [z_coords(5:9); z_coords(15)];
rho_disks      = [2830;2830;2870;4430;4430;8000];
thick_disks    = [0.02;0.02;0.02;0.02;0.02;0.03];
Dext_disks     = [0.19;0.21;0.23;0.24;0.25;0.25];
d_extShaft_dsk = [0.08;0.08;0.08;0.08;0.08;0.055];
m0             = [0;0;0;0;0;1e-4/0.125];

%% Bearing positions
z_bearings = [z_coords(2); z_coords(13)];

%% Mesh (Timoshenko type 3, element 10 refined)
zero_vec         = zeros(numEl,1);
n_intNodes       = zeros(numEl,1); n_intNodes(10) = 6;
isGraded         = zeros(numEl,1);
shaftElType      = 3*ones(numEl,1);
shaftElProperties = [rho, E, poisson];
d_shaftEl        = [d_IntShaft, d_ExtShaft, d_IntShaft, d_ExtShaft];
shaftLoad        = zeros(numEl,2);

rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, shaftElType, shaftElProperties, shaftLoad, beta);

%% Disk and bearing assembly
z_nodes_mesh = [rotor.nodes.node]; z_coord_mesh = [rotor.nodes.coord];

[found_flag, idx_disks] = ismember(z_disks, z_coord_mesh);
idx_disks = idx_disks(found_flag);

idx_b1 = find(z_coord_mesh == z_bearings(1), 1);
idx_b2 = find(z_coord_mesh == z_bearings(2), 1);

rotor.disk(6) = struct('type',[],'node',[],'thick',[],'D_ext',[],'D_int',[],'rho',[],'m0',[]);
[rotor.disk.type]  = deal(1);
[rotor.disk.node]  = num2cell(idx_disks){:};
[rotor.disk.thick] = num2cell(thick_disks){:};
[rotor.disk.D_ext] = num2cell(Dext_disks){:};
[rotor.disk.D_int] = num2cell(d_extShaft_dsk){:};
[rotor.disk.rho]   = num2cell(rho_disks){:};
[rotor.disk.m0]    = num2cell(m0){:};

% front bearing: constant stiffness, speed-dependent damping
rotorSpeed_pts = linspace(0,500,7)*2*pi;
kxx1 = 1.37e7*ones(1,7); cxx1 = [39.3 33.4 19.4 13.5 9.3 7.0 5.1]; z0 = zeros(1,7);
rotor.bearing(1) = struct('type',13,'node',z_nodes_mesh(idx_b1),...
    'Speed_points',rotorSpeed_pts,'kxx_points',kxx1,'kxy_points',z0,...
    'kyx_points',z0,'kyy_points',kxx1,'cxx_points',cxx1,'cxy_points',z0,...
    'cyx_points',z0,'cyy_points',cxx1);

% rear bearing: speed-dependent stiffness and damping
kxx2 = [2.94 2.94 5.30 7.16 8.44 10.5 12.5]*1e6;
cxx2 = [9719.7 9719.7 7855.7 6055.9 4216.1 3352.2 2441.2];
rotor.bearing(2) = struct('type',13,'node',z_nodes_mesh(idx_b2),...
    'Speed_points',rotorSpeed_pts,'kxx_points',kxx2,'kxy_points',z0,...
    'kyx_points',z0,'kyy_points',kxx2,'cxx_points',cxx2,'cxy_points',z0,...
    'cyx_points',z0,'cyy_points',cxx2);

rotor.forcing(1) = struct('type',1);

%% Rotor schematic
figureRotor(rotor);

%% Campbell diagram
rotorSpeed = linspace(0,500,10)*2*pi;
[natPuls, Mode, kappa] = charRoots(rotor, rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);
selectModes   = 1:6;
NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0   = kappa_sort(:,selectModes,:);
NX = [1 0.5]; isDamped = true;
plotCampbell(rotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort0);

%% Critical speeds and mode shapes
method = 2; NX_cs = 1;
[criticalSpeeds, modeCrits] = critSpeeds(rotor, NX_cs, isDamped, method, "num_crit", 6);

figure('Name','Critical Modes','NumberTitle','off')
for jj = 1:6
    subplot(2,3,jj)
    plotMode(rotor, modeCrits(:,jj), 'number', jj, ...
        'frequency', criticalSpeeds(jj)/(2*pi), ...
        'rotorSpeed', criticalSpeeds(jj)*30/pi);
    fprintf('Critical Speed %d: %.3f rpm\n', jj, criticalSpeeds(jj)*30/pi);
end

%% Stability: root locus and log decrement
plotRootLocus(rotorSpeed*30/pi, NatPuls_sort(2,:))
plotRootLocus(rotorSpeed*30/pi, NatPuls_sort(4,:))
deltaLog = LogDecrement(NatPuls_sort);
plotLogDecrement(rotorSpeed/(2*pi), NatPuls_sort(1,:))

%% Unbalance frequency response
rotorSpeed_frf = (0:0.1:500)*2*pi;
response = FRF(rotor, rotorSpeed_frf);
plotFRF(rotorSpeed_frf*30/pi, response, [5.1, 5.2])
plotFRF(rotorSpeed_frf*30/pi, response, [21.1, 21.2])

% deflected shape at 4th critical speed
[~, idx_crit] = min(abs(rotorSpeed_frf - criticalSpeeds(4)));
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed_frf(idx_crit)*30/pi)

%% Run-up simulation: 0 → 500 Hz in 10 s
alpha   = [0.5*100*pi, 0, 0];   % phi = 0.5*100*pi*t^2 → Omega_end = 100*pi = 500 Hz * 2pi
t_span  = [0, 10];

numEl_mesh = length(z_nodes_mesh);
initPos = zeros(4*numEl_mesh, 1);
initVel = zeros(4*numEl_mesh, 1);

[time, q, qdot, speed] = runUp(rotor, t_span, alpha, initPos, initVel, nr=12);

%% Post-processing: response at node 21
node_target = 21;
idx_u = (node_target-1)*4 + 1;
idx_v = (node_target-1)*4 + 2;

disp_v     = q(:, idx_v);
disp_total = sqrt(q(:,idx_u).^2 + q(:,idx_v).^2);
rpm_axis   = speed * 30 / pi;

% transient response vs time
figure('Name','Runup_Time','NumberTitle','off','Color','w');
plot(time, disp_v, 'Color',[0 0.4470 0.7410],'LineWidth',1.5)
grid on; grid minor; box on;
xlabel('Time t [s]'); ylabel('Vertical Displacement v [m]');
title(['Transient Response — Node ' num2str(node_target)]);
xlim([0, max(time)]);

% unbalance amplitude vs speed
figure('Name','Runup_RPM','NumberTitle','off','Color','w');
plot(rpm_axis, disp_total,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
grid on; grid minor; box on;
xlabel('Rotor Speed \Omega [rpm]'); ylabel('Total Displacement Amplitude r [m]');
title(['Unbalance Response (Run-up) — Node ' num2str(node_target)]);

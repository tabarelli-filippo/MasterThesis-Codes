% ES6_11MOD  Modal and frequency-domain analysis of a multi-disk rotor on
%   three fluid-film journal bearings (Exercise 6.11 — modal analysis).
%
%   ROTOR DESCRIPTION
%     Shaft : L = 2.8 m, 28 Timoshenko elements (type 2, Cowper),
%             variable diameter (d_shaft array, [mm])
%     Disks : 5 disks at z = 0.4, 0.8, 1.6, 2.0, 2.4 m
%             D_ext = 200 mm, D_int = 110 mm; thicknesses vary
%             Unbalance: 5 g on disk 4 (z = 2.0 m)
%     Bearings: 3 short fluid-film bearings (type 7) at z = 0, 1.2, 2.8 m
%               L = 30 mm, D = 100 mm, c = 2 mm, eta = 0.030 Pa·s
%               Static loads in ratio 1:1:1 (equal load sharing)
%
%   ANALYSIS PERFORMED
%     1. Total rotor mass and schematic (massRotor, figureRotor)
%     2. Critical speeds (iterative method 3 with initial estimates,
%        isDamped = true) — 6 critical speeds
%     3. Natural frequencies and Campbell diagram (sorted with MAC):
%        two plots (unsorted and sorted modes) for comparison
%     4. Mode shape comparison at two speeds (plotMode)
%     5. Orbit plots (plotOrbit)
%     6. Unbalance response at nodes 10 and 21 (FRF, plotFRF)
%     7. Mode shapes at critical speeds (plotMode, 2×3 layout)
%     8. Orbit plots at critical speeds
%     9. Root locus for one mode (plotRootLocus)
%    10. Stress analysis (stressAnalysis, plotStressAnalysis)
%
% SEE ALSO
%   es6_11ODE, es6_10mod, charRoots, critSpeeds, sortModesMAC, FRF

clear
close all
addpath(genpath('auxFunc/'))

%% Rotor parameters
L_shaft = 2.8;           % shaft length [m]
numel   = 28;            % number of elements
L_El    = L_shaft / numel;

z_bears  = [0, 1.2, 2.8];    % bearing positions [m]
node_bears = int8(z_bears / L_El + 1);

z_disks  = [0.4 0.8 1.6 2.0 2.4];  % disk positions [m]
node_disks  = int8(z_disks / L_El + 1);
Dext_disks  = 0.2;           % disk external diameter [m]
Dint_disks  = 0.11;          % disk internal diameter [m]
thick_disks = [25 25 25 100 25]*1e-3;  % disk thicknesses [m]
m0_disks = [0 0 0 5e-3 0];   % unbalance masses [kg]
m0_phase = [0 0 0 0 0];      % unbalance phases [rad]

rho     = 7800;          % density [kg/m³]
E       = 200e9;         % Young's modulus [Pa]
poisson = 0.285;
G       = E / 2 / (1 + poisson);

% shaft diameter at each node [m]
d_shaft = [100 100 38 110 38 110 38 110 38 38 38 100 100 38 38 38 110 38 38 ...
           110 110 38 38 38 110 38 100 100]*1e-3;

%% Rotor structure definition
rotor.nodes(numel+1) = struct('node',[],'coord',[]);
rotor.shaft(numel)   = struct('type',[],'node1',[],'node2',[],'d_int1',[],...
    'd_ext1',[],'d_int2',[],'d_ext2',[],'rho',[],'E',[],'G',[]);
rotor.disk(5)    = struct('type',[],'node',[],'thick',[],'D_int',[],...
    'D_ext',[],'rho',[],'m0',[],'gamma',[]);
rotor.bearing(3) = struct('type',[],'node',[]);

% nodes
nodes   = num2cell(1:numel+1);
z_nodes = num2cell(0:L_El:2.8);
rotor.nodes = struct('node', nodes, 'coord', z_nodes);

% disks
[rotor.disk.type]  = deal(1);
[rotor.disk.rho]   = deal(rho);
[rotor.disk.D_int] = deal(Dint_disks);
[rotor.disk.D_ext] = deal(Dext_disks);
[rotor.disk.node]  = num2cell(node_disks){:};
[rotor.disk.thick] = num2cell(thick_disks){:};
[rotor.disk.m0]    = num2cell(m0_disks){:};
[rotor.disk.gamma] = num2cell(m0_phase){:};

% shaft (Timoshenko type 2, variable diameter)
[rotor.shaft.type]  = deal(2);
[rotor.shaft.rho]   = deal(rho);
[rotor.shaft.E]     = deal(E);
[rotor.shaft.G]     = deal(G);
[rotor.shaft.d_int1] = deal(0);
[rotor.shaft.d_int2] = deal(0);
node1_cell = num2cell(1:numel); node2_cell = num2cell(2:numel+1);
[rotor.shaft.node1] = node1_cell{:};
[rotor.shaft.node2] = node2_cell{:};
Dext_shaft_cell = num2cell(d_shaft);
[rotor.shaft.d_ext1] = Dext_shaft_cell{:};
[rotor.shaft.d_ext2] = Dext_shaft_cell{:};

% fluid-film bearings (type 7), equal load sharing
disp('Fluid-film bearings — load ratio 1:1:1')
L_b = 0.03; D_b = 0.1; c_b = 2e-3*D_b; eta_b = 0.030;
W   = 136.0976 * 9.81;
F   = W / 3 * [1, 1, 1];
[rotor.bearing.type] = deal(7);
[rotor.bearing.node] = num2cell(node_bears){:};
[rotor.bearing.F]    = num2cell(F){:};
[rotor.bearing.L]    = deal(L_b);
[rotor.bearing.D]    = deal(D_b);
[rotor.bearing.c]    = deal(c_b);
[rotor.bearing.eta]  = deal(eta_b);

isDamped = true;
[rotor.forcing.type] = deal(1);

%% Mass and schematic
mass = massRotor(rotor);
figureRotor(rotor);

%% Critical speeds (iterative method 3 with initial estimates)
NX = 1; method = 3; num_crit = 6; max_iter = 40; toll = 1e-3;
initialSpeeds = [600;700;1000;1050;1800;2210;2350;2720]*pi/30;
[critical_speeds, mode_shapes] = critSpeeds(rotor, NX, isDamped, method, ...
    "num_crit", num_crit, "max_iter", max_iter, "toll", toll, ...
    "initial_estimates", initialSpeeds);
for kk = 1:6
    fprintf('Critical Speed %d: %.3f rpm\n', kk, critical_speeds(kk)*30/pi);
end

%% Campbell diagram (unsorted and MAC-sorted)
RotorSpeed = (10:100:3000)*pi/30;
[natPuls, Mode, kappa] = charRoots(rotor, RotorSpeed);
[natPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

NatPuls0       = natPuls(1:7,:);      kappa0       = kappa(:,1:7,:);
NatPuls_sort0  = natPuls_sort(1:7,:); kappa_sort0  = kappa_sort(:,1:7,:);

plotCampbell(RotorSpeed, NatPuls0, NX, isDamped, kappa0);        % unsorted
plotCampbell(RotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort0); % sorted

%% Mode shapes at two speeds
figure('Name','Modes','NumberTitle','off')
subplot(2,2,1); plotMode(rotor,Mode_sort(:,1,1),'number',1,'frequency',abs(natPuls_sort(1,1))/2/pi,'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,2); plotMode(rotor,Mode_sort(:,3,1),'number',3,'frequency',abs(natPuls_sort(3,1))/2/pi,'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,3); plotMode(rotor,Mode_sort(:,1,15),'number',1,'frequency',abs(natPuls_sort(1,15))/2/pi,'rotorSpeed',RotorSpeed(15)*30/pi);
subplot(2,2,4); plotMode(rotor,Mode_sort(:,3,15),'number',3,'frequency',abs(natPuls_sort(3,15))/2/pi,'rotorSpeed',RotorSpeed(15)*30/pi);

%% Orbit plots
figure('Name','Orbits','NumberTitle','off')
subplot(1,2,1); plotOrbit(Mode_sort(:,1,1), [1,5,10],  RotorSpeed(1)/2/pi, 'titletext','Orbit Nodes 1,5,10 — low speed');
subplot(1,2,2); plotOrbit(Mode_sort(:,1,15),[1,5,10],  RotorSpeed(15)/2/pi,'titletext','Orbit Nodes 1,5,10 — mid speed');

%% Unbalance frequency response
RotorSpeed = (10:10:3000)*pi/30;
response = FRF(rotor, RotorSpeed);
plotFRF(RotorSpeed*30/pi, response, [10.1 21.1])

%% Mode shapes at critical speeds (2x3 layout)
figure('Name','Critical Modes','NumberTitle','off')
for ii = 1:6
    subplot(2,3,ii)
    plotMode(rotor, mode_shapes(:,ii), 'number', ii, ...
        'frequency', critical_speeds(ii)/(2*pi), ...
        'rotorSpeed', critical_speeds(ii)*30/pi);
    title(sprintf('Mode %d — %.0f rpm', ii, critical_speeds(ii)*30/pi));
end

%% Orbit plots at critical speeds
figure('Name','Critical Orbits','NumberTitle','off')
subplot(1,2,1); plotOrbit(mode_shapes(:,1),[1,5,10],critical_speeds(1)/2/pi,'titletext','Mode 1 orbit');
subplot(1,2,2); plotOrbit(mode_shapes(:,3),[1,5,10],critical_speeds(3)/2/pi,'titletext','Mode 3 orbit');

%% Root locus (first mode)
RotorSpeed_rl = 10:10:3000;
plotRootLocus(RotorSpeed_rl*30/pi, NatPuls_sort0(1,:));

%% Stress analysis
[nu_static, ~, ~] = stressAnalysis(rotor, response(:,100), 350e6);
plotStressAnalysis(rotor, response(:,100), nu_static);

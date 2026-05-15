% ES6_10MOD  Modal analysis of a simple rotor with three disks and two
%   pinned bearings (Exercise 6.10 — modal/frequency-domain analysis).
%
%   ROTOR DESCRIPTION
%     Shaft : L = 1.6 m, d = 75 mm, solid, steel (E = 200 GPa, rho = 7800 kg/m³)
%     Disks : three identical disks (D_ext = 400 mm, thick = 80 mm)
%             at z = 0.0, 0.8, 1.2 m; each with a 20 g unbalance mass
%     Bearings: two pinned supports (type 1) at z = 0.4 m and z = 1.6 m
%     Mesh  : 8 Euler-Bernoulli elements (type 1)
%
%   ANALYSIS PERFORMED
%     1. Rotor schematic (figureRotor)
%     2. Natural frequencies vs speed (charRoots) and Campbell diagram
%        with 1X, 2X, 3X engine-order lines (plotCampbell)
%     3. Critical speeds for 1X, 2X, 3X excitation using both direct
%        (method 1) and iterative (method 2) solvers (critSpeeds)
%     4. Unbalance frequency response at nodes 1 and 5 (FRF, plotFRF)
%     5. Mode shape visualisation at Omega = 0 and Omega = 100*pi/30 rad/s
%        (plotMode)
%     6. Orbit plots at two node pairs (plotOrbit)
%     7. Mode shapes at critical speeds (plotMode)
%     8. Root locus for the first mode (plotRootLocus)
%     9. Stress analysis at one speed step (stressAnalysis, plotStressAnalysis)
%
% SEE ALSO
%   es6_10ODE, es6_11mod, charRoots, critSpeeds, FRF, plotCampbell,
%   plotMode, plotOrbit, plotRootLocus, stressAnalysis

clear
close all
addpath(genpath('auxFunc/'))

%% Rotor parameters
L_shaft   = 1.6;       % shaft length [m]
d_shaft   = 75e-3;     % shaft diameter [m]
D_disk    = 0.4;       % disk external diameter [m]
thick_disk = 0.08;     % disk thickness [m]
z_disks   = [0.0 0.8 1.2];   % disk axial positions [m]
z_bearings = [0.4 1.6];      % bearing axial positions [m]
E   = 200e9;           % Young's modulus [Pa]
v_poisson = 0.27;
G   = E / 2 / (1 + v_poisson);
rho = 7800;            % density [kg/m³]
Rotor_speed = 3000 * pi / 30;  % operating speed [rad/s]

numEl = 8;             % number of shaft finite elements
L_El  = L_shaft / numEl;

%% Rotor structure definition
rotor.nodes(numEl+1) = struct('node',[],'coord',[]);
rotor.shaft(numEl)   = struct('type',[],'node1',[],'node2',[],...
    'd_int1',[],'d_ext1',[],'d_int2',[],'d_ext2',[],'rho',[],'E',[],'G',[]);
rotor.disk(3)    = struct('type',[],'node',[],'thick',[],'D_int',[],...
    'D_ext',[],'rho',[],'m0',[],'gamma',[]);
rotor.bearing(2) = struct('type',[],'node',[]);
rotor.forcing(1) = struct('type',[]);

% nodes
nodes   = num2cell(1:numEl+1);
z_nodes = num2cell(0:L_El:L_shaft);
rotor.nodes = struct('node', nodes, 'coord', z_nodes);

% shaft elements (Euler-Bernoulli, solid, uniform)
[rotor.shaft.type]  = deal(1);
[rotor.shaft.rho]   = deal(rho);
[rotor.shaft.E]     = deal(E);
[rotor.shaft.G]     = deal(G);
[rotor.shaft.d_int1] = deal(0);
[rotor.shaft.d_int2] = deal(0);
[rotor.shaft.d_ext1] = deal(d_shaft);
[rotor.shaft.d_ext2] = deal(d_shaft);
node1_cell = num2cell(1:numEl);   node2_cell = num2cell(2:numEl+1);
[rotor.shaft.node1] = node1_cell{:};
[rotor.shaft.node2] = node2_cell{:};

% disks
[rotor.disk.type]  = deal(1);
[rotor.disk.rho]   = deal(rho);
[rotor.disk.D_int] = deal(d_shaft);
[rotor.disk.D_ext] = deal(D_disk);
[rotor.disk.thick] = deal(thick_disk);
node_disks      = int8(z_disks / L_El + 1);
[rotor.disk.node]  = num2cell(node_disks){:};
% unbalance: auxiliary mass m0 on each disk
m0_disks  = [20e-3 20e-3 20e-3]; % [kg]
m0_phase  = [0 0 0];             % [rad]
[rotor.disk.m0]    = num2cell(m0_disks){:};
[rotor.disk.gamma] = num2cell(m0_phase){:};

% bearings (pinned, type 1)
node_bearings = int8(z_bearings/L_El + 1);
[rotor.bearing.type] = deal(1);
[rotor.bearing.node] = num2cell(node_bearings){:};

% forcing: mass unbalance
[rotor.forcing.type] = deal(1);

%% Rotor schematic
figureRotor(rotor);

%% Natural frequencies and Campbell diagram
RotorSpeed = (0:100:5000)*pi/30;  % speed sweep [rad/s]
[natPuls, Mode, kappa] = charRoots(rotor, RotorSpeed);
natPuls0 = natPuls(1:5,:);

NX       = [1 2 3];
isDamped = false;
kappa0   = kappa(:,1:2:10,:);
plotCampbell(RotorSpeed, natPuls0, NX, isDamped, kappa0)

%% Critical speeds (both methods, 1X–3X excitation)
isDamped  = false;
num_crit  = 5;
max_iter  = 20;
toll      = 1e-3;

RotorSpeed = (0:10:5000)*pi/30;
methods    = [1, 2];
labels     = {"Direct", "Iterative"};

for m = 1:length(methods)
    fprintf('\n=============================================');
    fprintf('\nMETHOD: %s', labels{m});
    fprintf('\n=============================================\n');
    for N_exc = 1:3
        [crit_speeds, ~] = critSpeeds(rotor, N_exc, isDamped, methods(m), ...
            "num_crit", num_crit, "max_iter", max_iter, "toll", toll);
        fprintf('\n--- %dX Excitation ---', N_exc);
        for jj = 1:num_crit
            fprintf('\n%dX - Critical Speed %d: %.3f [rpm]', N_exc, jj, crit_speeds(jj)*30/pi);
        end
        fprintf('\n');
    end
end

%% Unbalance frequency response
response = FRF(rotor, RotorSpeed);
plotFRF(RotorSpeed*30/pi, response, [1.1 5.1])

%% Mode shape comparison (two speeds)
figure('Name','Modes','NumberTitle', 'off')
subplot(2,2,1)
plotMode(rotor, Mode(:,1,1),  'number',1,'frequency',abs(natPuls0(1,1))/2/pi, 'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,2)
plotMode(rotor, Mode(:,3,1),  'number',2,'frequency',abs(natPuls0(3,1))/2/pi, 'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,3)
plotMode(rotor, Mode(:,1,25), 'number',3,'frequency',abs(natPuls0(1,25))/2/pi,'rotorSpeed',RotorSpeed(25)*30/pi);
subplot(2,2,4)
plotMode(rotor, Mode(:,3,25), 'number',4,'frequency',abs(natPuls0(3,25))/2/pi,'rotorSpeed',RotorSpeed(25)*30/pi);

%% Orbit plots
figure('Name','Orbits','NumberTitle', 'off')
subplot(1,2,1)
plotOrbit(Mode(:,1,25), [1,2], natPuls0(1,2), 'titletext','Orbit Nodes 1 and 2');
subplot(1,2,2)
plotOrbit(Mode(:,1,25), [1,5], natPuls0(1,2), 'titletext','Orbit Nodes 1 and 5');

%% Mode shapes at critical speeds
% (requires mode_shapes1 and critical_speeds1 from critSpeeds with nargout=2)
figure('Name','Critical Modes','NumberTitle', 'off')
subplot(2,2,1)
plotMode(rotor, mode_shapes1(:,1),'number',1,'frequency',abs(critical_speeds1(1))/2/pi,'rotorSpeed',critical_speeds1(1)*30/pi);
subplot(2,2,2)
plotMode(rotor, mode_shapes1(:,2),'number',2,'frequency',abs(critical_speeds1(2))/2/pi,'rotorSpeed',critical_speeds1(2)*30/pi);
subplot(2,2,3)
plotMode(rotor, mode_shapes1(:,3),'number',3,'frequency',abs(critical_speeds1(3))/2/pi,'rotorSpeed',critical_speeds1(3)*30/pi);
subplot(2,2,4)
plotMode(rotor, mode_shapes1(:,4),'number',4,'frequency',abs(critical_speeds1(4))/2/pi,'rotorSpeed',critical_speeds1(4)*30/pi);

%% Root locus (first mode)
RotorSpeed = (0:100:5000)*pi/30;
plotRootLocus(RotorSpeed*30/pi, natPuls0(1,:));

%% Stress analysis at one speed
[nu_static, ~, ~] = stressAnalysis(rotor, response(:,400), 350e6);
plotStressAnalysis(rotor, response(:,400), nu_static);

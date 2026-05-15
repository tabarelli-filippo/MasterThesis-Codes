% ES6_10ODE  Nonlinear time-domain simulation of a rotor with three disks
%   and two pinned bearings (Exercise 6.10 — time-integration analysis).
%
%   Uses the same rotor geometry as es6_10mod but with Timoshenko beam
%   elements (type 3) and a single unbalance mass on disk 1 only. The
%   steady-state orbit is analysed via Poincaré map and bifurcation plot.
%
%   ROTOR DESCRIPTION
%     Shaft : L = 1.6 m, d = 75 mm, solid, steel
%     Disks : three identical disks at z = 0.0, 0.8, 1.2 m
%             Unbalance: 20 g on disk 1 only (z = 0.0 m)
%     Bearings: two pinned supports (type 1) at z = 0.4 m and z = 1.6 m
%     Mesh  : 8 Timoshenko elements (type 3, Hutchinson shear coefficient)
%
%   ANALYSIS PERFORMED
%     1. Rotor schematic (figureRotor)
%     2. Time simulation from t = 0 to t = 400 s at Omega = 3000 rpm
%        with 4-mode reduction (timeSimulation, nr = 4)
%     3. Phase-plane orbit plot for nodes 1–2 (steady-state tail)
%     4. Poincaré map at the operating speed (poincareMap)
%     5. Bifurcation diagram vs speed (plotBifurcation)
%
%   NOTE: The long simulation time (t = 400 s) is chosen to ensure the
%   transient has fully decayed and only the steady-state attractor remains.
%   The initial condition uses a small non-zero x-displacement (1e-5 m) at
%   all nodes to break symmetry.
%
% SEE ALSO
%   es6_10mod, timeSimulation, poincareMap, plotBifurcation

clear
close all
addpath("auxFunc/");

%% Rotor parameters
L_shaft   = 1.6;
d_shaft   = 75e-3;
D_disk    = 0.4;
thick_disk = 0.08;
z_disks   = [0.0 0.8 1.2];
z_bearings = [0.4 1.6];
E   = 200e9;
v_poisson = 0.27;
G   = E / 2 / (1 + v_poisson);
rho = 7800;
Rotor_speed = 3000 * pi / 30;  % [rad/s]

numEl = 8;
L_El  = L_shaft / numEl;

%% Rotor structure definition
rotor.nodes(numEl+1) = struct('node',[],'coord',[]);
rotor.shaft(numEl)   = struct('type',[],'node1',[],'node2',[],...
    'd_int1',[],'d_ext1',[],'d_int2',[],'d_ext2',[],'rho',[],'E',[],'G',[]);
rotor.disk(3)    = struct('type',[],'node',[],'thick',[],'D_int',[],...
    'D_ext',[],'rho',[],'m0',[],'gamma',[]);
rotor.bearing(2) = struct('type',[],'node',[]);
rotor.forcing(1) = struct('type',[]);

nodes   = num2cell(1:numEl+1);
z_nodes = num2cell(0:L_El:L_shaft);
rotor.nodes = struct('node', nodes, 'coord', z_nodes);

% shaft elements (Timoshenko type 3 — Hutchinson shear coefficient)
[rotor.shaft.type]   = deal(3);
[rotor.shaft.rho]    = deal(rho);
[rotor.shaft.E]      = deal(E);
[rotor.shaft.G]      = deal(G);
[rotor.shaft.d_int1] = deal(0);
[rotor.shaft.d_int2] = deal(0);
[rotor.shaft.d_ext1] = deal(d_shaft);
[rotor.shaft.d_ext2] = deal(d_shaft);
node1_cell = num2cell(1:numEl); node2_cell = num2cell(2:numEl+1);
[rotor.shaft.node1] = node1_cell{:};
[rotor.shaft.node2] = node2_cell{:};

% disks
[rotor.disk.type]    = deal(1);
[rotor.disk.rho]     = deal(rho);
[rotor.disk.D_int]   = deal(d_shaft);
[rotor.disk.D_ext]   = deal(D_disk);
[rotor.disk.thick]   = deal(thick_disk);
node_disks = int8(z_disks / L_El + 1);
[rotor.disk.node]    = num2cell(node_disks){:};
% unbalance on disk 1 only
m0_disks = [20e-3 0 0];
m0_phase = [0 0 0];
[rotor.disk.m0]    = num2cell(m0_disks){:};
[rotor.disk.gamma] = num2cell(m0_phase){:};

% bearings (pinned)
node_bearings = int8(z_bearings/L_El + 1);
[rotor.bearing.type] = deal(1);
[rotor.bearing.node] = num2cell(node_bearings){:};

% forcing: mass unbalance
[rotor.forcing.type] = deal(1);

%% Rotor schematic
figureRotor(rotor);

%% Time-domain simulation
t_span  = [0, 400];              % integration interval [s]
initPos = zeros(4*(numEl+1), 1);
initPos(1:4:end) = 1e-5;         % small initial x-displacement to break symmetry
initVel = zeros(4*(numEl+1), 1);

[time, q, qdot] = timeSimulation(rotor, t_span, Rotor_speed, initPos, initVel, nr=4);

%% Post-processing: steady-state tail (t > 390 s)
idx = time > 390;

% phase-plane orbit (node 1: DOFs 1 and 2)
figure('Name','Steady-state Orbit','NumberTitle','off');
plot(q(idx,1), q(idx,2));
xlabel('u [m]'); ylabel('v [m]');
title('Steady-state orbit — Node 1 (x vs y)');
grid on; axis equal;

% Poincaré map
poincareMap(time(idx), q(idx,1), q(idx,2), Rotor_speed)

% bifurcation diagram
plotBifurcation(Rotor_speed, time(idx), q(idx,1))

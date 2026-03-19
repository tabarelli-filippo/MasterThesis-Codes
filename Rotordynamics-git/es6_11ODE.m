% ES6_11ODE  Nonlinear time-domain simulation of a multi-disk rotor on
%   three nonlinear fluid-film bearings (Exercise 6.11 — time integration).
%
%   Uses the same rotor geometry as es6_11mod but with nonlinear bearing
%   type 7.1 (Ocvirk short-bearing, pi-film model) and an unbalance on two
%   disks. The simulation runs at Omega = 1000 rpm and produces the
%   steady-state orbit at a selected node.
%
%   ROTOR DESCRIPTION
%     Shaft : L = 2.8 m, 28 Timoshenko elements (type 2)
%     Disks : 5 disks; unbalances on disks 1 and 4 (5 g each)
%     Bearings: 3 nonlinear fluid-film bearings (type 7.1) at z = 0, 1.2, 2.8 m
%               L = 30 mm, D = 110 mm, c = 0.5*D, eta = 0.030 Pa·s
%               Static loads proportional to rotor weight (equal sharing)
%
%   ANALYSIS PERFORMED
%     1. Rotor schematic (figureRotor)
%     2. Time simulation t = [0, 51] s at Omega = 1000 rpm
%        with 5-mode reduction (timeSimulation, nr = 5)
%     3. Phase-plane orbit at node 10 (DOFs 40–41)
%
%   NOTE: The large clearance c = 0.5*D is intentionally exaggerated to
%   demonstrate nonlinear behaviour within a short simulation window.
%
% SEE ALSO
%   es6_11mod, timeSimulation, nonLinBearingMatrix, poincareMap

clear
close all
addpath(genpath('auxFunc/'));

%% Rotor parameters
L_shaft = 2.8; numEl = 28; L_El = L_shaft / numEl;
z_bears = [0, 1.2, 2.8]; node_bears = int8(z_bears / L_El + 1);
z_disks = [0.4 0.8 1.6 2.0 2.4];
node_disks = int8(z_disks / L_El + 1);
Dext_disks = 0.2; Dint_disks = 0.11;
thick_disks = [25 25 25 100 25]*1e-3;
m0_disks = [5e-3 0 0 5e-3 0]; m0_phase = [0 0 0 0 0];
rho = 7800; E = 200e9; poisson = 0.285; G = E/2/(1+poisson);
d_shaft = [100 100 38 110 38 110 38 110 38 38 38 100 100 38 38 38 110 38 38 ...
           110 110 38 38 38 110 38 100 100]*1e-3;
RotorSpeed = 1000*pi/30;  % operating speed [rad/s]

%% Rotor structure definition
rotor.nodes(numEl+1) = struct('node',[],'coord',[]);
rotor.shaft(numEl)   = struct('type',[],'node1',[],'node2',[],'d_int1',[],...
    'd_ext1',[],'d_int2',[],'d_ext2',[],'rho',[],'E',[],'G',[]);
rotor.disk(5)    = struct('type',[],'node',[],'thick',[],'D_int',[],...
    'D_ext',[],'rho',[],'m0',[],'gamma',[]);
rotor.bearing(3) = struct('type',[],'node',[]);

nodes = num2cell(1:numEl+1); z_nodes = num2cell(0:L_El:2.8);
rotor.nodes = struct('node', nodes, 'coord', z_nodes);

[rotor.disk.type] = deal(1); [rotor.disk.rho] = deal(rho);
[rotor.disk.D_int] = deal(Dint_disks); [rotor.disk.D_ext] = deal(Dext_disks);
[rotor.disk.node]  = num2cell(node_disks){:};
[rotor.disk.thick] = num2cell(thick_disks){:};
[rotor.disk.m0]    = num2cell(m0_disks){:};
[rotor.disk.gamma] = num2cell(m0_phase){:};

[rotor.shaft.type]   = deal(2);
[rotor.shaft.rho]    = deal(rho); [rotor.shaft.E] = deal(E); [rotor.shaft.G] = deal(G);
[rotor.shaft.d_int1] = deal(0);  [rotor.shaft.d_int2] = deal(0);
node1_cell = num2cell(1:numEl); node2_cell = num2cell(2:numEl+1);
[rotor.shaft.node1] = node1_cell{:}; [rotor.shaft.node2] = node2_cell{:};
Dext_shaft_cell = num2cell(d_shaft);
[rotor.shaft.d_ext1] = Dext_shaft_cell{:};
[rotor.shaft.d_ext2] = Dext_shaft_cell{:};

% nonlinear fluid-film bearings (type 7.1)
L_b = 0.03; D_b = 0.11; c_b = 0.5*D_b; eta_b = 0.03;
mass_rotor = massRotor(rotor);
F = mass_rotor*9.81/3 * [1,1,1];
[rotor.bearing.type] = deal(7.1);
[rotor.bearing.node] = num2cell(node_bears){:};
[rotor.bearing.F]    = num2cell(F){:};
[rotor.bearing.L]    = deal(L_b); [rotor.bearing.D] = deal(D_b);
[rotor.bearing.c]    = deal(c_b); [rotor.bearing.eta] = deal(eta_b);
[rotor.forcing.type] = deal(1);

%% Rotor schematic
figureRotor(rotor);

%% Time-domain simulation
t_span  = [0, 51];
initPos = zeros(4*(numEl+1), 1);
initPos(2:4:end) = -1e-5;     % small initial y-displacement
initVel = zeros(4*(numEl+1), 1);

[time, q, qdot] = timeSimulation(rotor, t_span, RotorSpeed, initPos, initVel, nr=5);

%% Phase-plane orbit at node 10 (DOFs 40–41)
figure('Name','Orbit Node 10','NumberTitle','off');
plot(q(:,40), q(:,41));
xlabel('u [m]'); ylabel('v [m]');
title(sprintf('Orbit — Node 10, Omega = %.0f rpm', RotorSpeed*30/pi));
grid on; axis equal;

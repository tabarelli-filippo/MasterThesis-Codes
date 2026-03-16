clear
close all

addpath(genpath('auxFunc/'))
%% esercizio 6.10
% dati
L_shaft = 1.6; %[m]
d_shaft = 75e-3; %[m]

D_disk = 0.4; % [m]
thick_disk = 0.08; % [m]
z_disks = [0.0 0.8 1.2]; %[m]

z_bearings = [0.4 1.6]; %[m]

E = 200e9; %[Pa]
v_poisson = 0.27; %[-]
G = E / 2 / (1 + v_poisson);
rho = 7800; % [Kg/m^3]

Rotor_speed = 3000 * pi / 30; %[rad/s]

% schematizzare albero con 8 elementi
Nat_Freq_Hz = zeros(5,1);
numEl = 8;
L_El = L_shaft / numEl;

%% definizione della struttura
% inizializzazione elementi
rotor.nodes(numEl+1) = struct('node',[],'coord',[]);
rotor.shaft(numEl) = struct('type', [], 'node1', [], 'node2', [], 'd_int1', ...
    [], 'd_ext1', [],'d_int2',[],'d_ext2',[],'rho', [], 'E', [], 'G', []);
rotor.disk(3) = struct('type', [], 'node', [], 'thick', [], 'D_int', ...
    [], 'D_ext', [], 'rho',[],'m0',[],'gamma',[]);
rotor.bearing(2) = struct('type', [], 'node', []);
rotor.forcing(1) = struct('type', []);


% definizione nodo
nodes = num2cell(1:1:(numEl+1));
z_nodes = num2cell(0:L_El:L_shaft);
rotor.nodes = struct('node', nodes, 'coord', z_nodes);

% definizione elementi albero
%elemento albero Eulero-Bernoulli pieno
[rotor.shaft.type] = deal(1); %elemento Eulero
[rotor.shaft.rho] = deal(rho);
[rotor.shaft.E] = deal(E);
[rotor.shaft.G] = deal(G);
[rotor.shaft.d_int1] = deal(0); % elemento pieno
[rotor.shaft.d_int2] = deal(0);
[rotor.shaft.d_ext1] = deal(d_shaft);
[rotor.shaft.d_ext2] = deal(d_shaft);
%nodi
node1_cell = num2cell(1:numEl);
node2_cell = num2cell(2:numEl+1);
[rotor.shaft.node1] = node1_cell{:};
[rotor.shaft.node2] = node2_cell{:};

%dischi
%costanti
[rotor.disk.type] = deal(1);
[rotor.disk.rho] = deal(rho);
[rotor.disk.D_int] = deal(d_shaft);
[rotor.disk.D_ext] = deal(D_disk);
[rotor.disk.thick] = deal(thick_disk);
%nodi
node_disks = int8(z_disks / L_El + 1);
node_disks_cell = num2cell(node_disks);
[rotor.disk.node] = node_disks_cell{:};
%sbilanciamento
%espresso come masse ausiliarie m0
m0_disks = [20e-3 20e-3 20e-3]; %[Kg]
m0_phase = [0 0 0]; %[rad]
m0_disks_cell = num2cell(m0_disks);
[rotor.disk.m0] = m0_disks_cell{:};
m0_phase_cell = num2cell(m0_phase);
[rotor.disk.gamma] = m0_phase_cell{:};

% cuscinetti auto-allineanti, non ho caratteristiche cuscinetti
% cuscinetti tipo 1 posizionati a z = 0.4 e z = 1.6
node_bearings = int8(z_bearings/L_El + 1);
node_bears_cell = num2cell(node_bearings);
[rotor.bearing.type] = deal(1);
[rotor.bearing.node] = node_bears_cell{:};

% sbilanciamento
[rotor.forcing.type] = deal(1);
%% rotor picture
figureRotor(rotor);
%% calcolo delle frequenze naturali
RotorSpeed = (0:100:5000)*pi/30;

[natPuls,Mode,kappa] = charRoots(rotor,RotorSpeed);
natPuls0 = natPuls(1:5,:);

% plot Campbell Diagram
NX = [1 2 3];
isDamped = 0;
kappa0 = kappa(:,1:2:10,:);

plotCampbell(RotorSpeed,natPuls0,NX,isDamped,kappa0)

%% calcolo delle frequenze critiche
N_unbalance = 1;
isDamped = 0;
num_crit = 5;
max_iter = 20;
toll = 1e-3;

RotorSpeed = (0:10:5000)*pi/30;
methods = [1, 2];
labels = {"Direct", " Iterative"};

for m = 1:length(methods)
    current_method = methods(m);
    fprintf('\n=============================================');
    fprintf('\nMETHOD: %s', labels{m});
    fprintf('\n=============================================\n');
    
    for N_exc = 1:3 
        [crit_speeds, ~] = critSpeeds(rotor, N_exc, isDamped, current_method, ...
            "num_crit", num_crit, "max_iter", max_iter, "toll", toll);
        fprintf('\n--- %dX Excitation ---', N_exc);
        for jj = 1:num_crit
            rpm_val = crit_speeds(jj) * 30/pi;
            fprintf('\n%dX - Critical Speed %d: %.3f [rpm]', N_exc, jj, rpm_val);
        end
        fprintf('\n');
    end
end
%% risposta allo sbilanciamento
response = FRF(rotor,RotorSpeed);
plotFRF(RotorSpeed*30/pi,response,[1.1 5.1])

%% mode comparison
figure('Name','Modes','NumberTitle', 'off')
subplot(2,2,1)
plotMode(rotor,Mode(:,1,1),'number',1,'frequency',abs(natPuls0(1,1))/2/pi,'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,2)
plotMode(rotor,Mode(:,3,1),'number',2,'frequency',abs(natPuls0(3,1))/2/pi,'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,3)
plotMode(rotor,Mode(:,1,25),'number',3,'frequency',abs(natPuls0(1,25))/2/pi,'rotorSpeed',RotorSpeed(25)*30/pi);
subplot(2,2,4)
plotMode(rotor,Mode(:,3,25),'number',4,'frequency',abs(natPuls0(3,25))/2/pi,'rotorSpeed',RotorSpeed(25)*30/pi);

%% rappresentazione delle orbite
figure('Name','Orbits','NumberTitle', 'off')
subplot(1,2,1)
plotOrbit(Mode(:,1,25),[1,2],natPuls0(1,2),'titletext','Orbit Nodes 1 and 2');
subplot(1,2,2)
plotOrbit(Mode(:,1,25),[1,5],natPuls0(1,2),'titletext','Orbit Node 1 ans 5');

%% Critical Modes
figure('Name','Modes','NumberTitle', 'off')
subplot(2,2,1)
plotMode(rotor,mode_shapes1(:,1),'number',1,'frequency',abs(critical_speeds1(1))/2/pi,'rotorSpeed',critical_speeds1(1)*30/pi);
subplot(2,2,2)
plotMode(rotor,mode_shapes1(:,2),'number',2,'frequency',abs(critical_speeds1(2))/2/pi,'rotorSpeed',critical_speeds1(2)*30/pi);
subplot(2,2,3)
plotMode(rotor,mode_shapes1(:,3),'number',3,'frequency',abs(critical_speeds1(3))/2/pi,'rotorSpeed',critical_speeds1(3)*30/pi);
subplot(2,2,4)
plotMode(rotor,mode_shapes1(:,4),'number',4,'frequency',abs(critical_speeds1(4))/2/pi,'rotorSpeed',critical_speeds1(4)*30/pi);

%% Root Locus: Pulsazione - 1
RotorSpeed = (0:100:5000)*pi/30;
plotRootLocus(RotorSpeed*30/pi,natPuls0(1,:));

%% Stress Analysis
[nu_static,nu_fatigue,N_life] = stressAnalysis(rotor,response(:,400),350e6);
plotStressAnalysis(rotor,response(:,400),nu_static);


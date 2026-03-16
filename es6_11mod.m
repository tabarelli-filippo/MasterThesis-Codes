clear
close all

addpath(genpath('auxFunc/'))
%% esercizio 6.11
% dati
L_shaft = 2.8; %[m]
numel = 28;
L_El = L_shaft / numel;

% posizione cuscinetti da sinistra
z_bears = [0, 1.2, 2.8]; %[m]
node_bears = int8(z_bears / L_El + 1);

% posizione dischi
z_disks = [0.4 0.8 1.6 2.0 2.4]; %[m]
node_disks = int8(z_disks / L_El + 1);
Dext_disks = 0.2; %[m]
Dint_disks = 0.11; %[m]
thick_disks = [25 25 25 100 25] * 1e-3; %[m]
m0_disks = [0 0 0 5e-3 0]; %[Kg]
m0_phase = [0 0 0 0 0]; %[rad]

rho = 7800; %[Kg/m^3]
E = 200e9; %[Pa]
poisson = 0.285; %[-]
G = E / 2 / (1 + poisson); %[Pa]

% shaft dimensions at each node
d_shaft = [100 100 38 110 38 110 38 110 38 38 38 100 100 38 38 38 110 38 38 ...
    110 110 38 38 38 110 38 100 100] * 1e-3; %[m]

RotorSpeed = 0 * pi / 30; %[rad/s]
%% caso 1 - cuscinetti isotropi
% inizializzazione strutture
rotor.nodes(numel+1) = struct('node',[],'coord',[]);
rotor.shaft(numel) = struct('type', [], 'node1', [], 'node2', [], 'd_int1', ...
    [], 'd_ext1', [],'d_int2',[],'d_ext2',[],'rho', [], 'E', [], 'G', []);
rotor.disk(5) = struct('type', [], 'node', [], 'thick', [], 'D_int', ...
    [], 'D_ext', [], 'rho',[],'m0',[],'gamma',[]);
rotor.bearing(3) = struct('type', [], 'node', []);


disp('Caso Cuscinetti a film fluido rapporto forze 1:1:1 - Rotor speed = 3000 [rpm]')
L = 0.03; %[m]
D = 0.1; %[m]
c = 2e-3 * D; %[m]
eta = 0.030; %[Pa s] viscosità olio
W = 136.0976 * 9.81; % peso del rotore
F = W / 3 * [1,1,1];
F_cell = num2cell(F);

% definizione del rotore:
%nodi
nodes = num2cell(1:1:(numel+1));
z_nodes = num2cell(0:L_El:2.8);
rotor.nodes = struct('node', nodes, 'coord', z_nodes);

%dischi
%costanti
[rotor.disk.type] = deal(1);
[rotor.disk.rho] = deal(rho);
[rotor.disk.D_int] = deal(Dint_disks);
[rotor.disk.D_ext] = deal(Dext_disks);
%nodi
node_disks_cell = num2cell(node_disks);
[rotor.disk.node] = node_disks_cell{:};
%spessori
thick_disks_cell = num2cell(thick_disks);
[rotor.disk.thick] = thick_disks_cell{:};
% sbilanciamento
m0_disks_cell = num2cell(m0_disks);
[rotor.disk.m0] = m0_disks_cell{:};
m0_phase_cell = num2cell(m0_phase);
[rotor.disk.gamma] = m0_phase_cell{:};

%albero
%costanti
[rotor.shaft.type] = deal(2); %elemento Timoshenko
[rotor.shaft.rho] = deal(rho);
[rotor.shaft.E] = deal(E);
[rotor.shaft.G] = deal(G);
[rotor.shaft.d_int1] = deal(0); % elemento pieno
[rotor.shaft.d_int2] = deal(0);
%nodi
node1_cell = num2cell(1:numel);
node2_cell = num2cell(2:numel+1);
[rotor.shaft.node1] = node1_cell{:};
[rotor.shaft.node2] = node2_cell{:};
%diametro ext
Dext_shaft_cell = num2cell(d_shaft);
[rotor.shaft.d_ext1] = Dext_shaft_cell{:};
[rotor.shaft.d_ext2] = Dext_shaft_cell{:};

% cuscinetti
[rotor.bearing.type] = deal(7); %cuscinetti tipo 7
node_bears_cell = num2cell(node_bears);
[rotor.bearing.node] = node_bears_cell{:};

%proprietà dei cuscinetti
[rotor.bearing.F] = F_cell{:};
[rotor.bearing.L] = deal(L);
[rotor.bearing.D] = deal(D);
[rotor.bearing.c] = deal(c);
[rotor.bearing.eta] = deal(eta);

isDamped = 1;

%forzante
[rotor.forcing.type] = deal(1);
%%
mass = massRotor(rotor);
figureRotor(rotor);
%% calcolo delle velocità critiche
NX = 1;
method = 3; num_crit = 6; max_iter = 40; toll=1e-3;

initialSpeeds = [600; 700; 1000; 1050; 1800; 2210; 2350; 2720];
initialSpeeds = initialSpeeds * pi/30;

[critical_speeds, mode_shapes] = critSpeeds(rotor, NX, isDamped, method, ...
    "num_crit", num_crit,"max_iter", max_iter, "toll", toll,"initial_estimates",initialSpeeds);

for kk = 1:6
    fprintf('Critical Speed %d : %.3f [rpm]\n',kk,critical_speeds(kk)*30/pi);
end

%% calcolo della risposta allo sbilanciamento
RotorSpeed = (10:100:3000)*pi/30; %[rad/s]
[natPuls,Mode,kappa] = charRoots(rotor,RotorSpeed);
[natPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

% non sorted
NatPuls0 = natPuls(1:1:7,:);
kappa0 = kappa(:,1:1:7,:);
isDamped = 1;
plotCampbell(RotorSpeed,NatPuls0,NX,isDamped,kappa0);

%sorted
NatPuls_sort0 = natPuls_sort(1:1:7,:);
kappa_sort0 = kappa_sort(:,1:1:7,:);
isDamped = 1;
plotCampbell(RotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);
%% modes comparison
figure('Name','Modes','NumberTitle', 'off')
subplot(2,2,1)
plotMode(rotor,Mode_sort(:,1,1),'number',1,'frequency',abs(natPuls_sort(1,1))/2/pi,'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,2)
plotMode(rotor,Mode_sort(:,3,1),'number',3,'frequency',abs(natPuls_sort(3,1))/2/pi,'rotorSpeed',RotorSpeed(1)*30/pi);
subplot(2,2,3)
plotMode(rotor,Mode_sort(:,1,15),'number',1,'frequency',abs(natPuls_sort(1,15))/2/pi,'rotorSpeed',RotorSpeed(15)*30/pi);
subplot(2,2,4)
plotMode(rotor,Mode_sort(:,3,15),'number',3,'frequency',abs(natPuls_sort(3,15))/2/pi,'rotorSpeed',RotorSpeed(15)*30/pi);

%% orbits comparisons
figure('Name','Orbits','NumberTitle', 'off')
subplot(1,2,1)
plotOrbit(Mode_sort(:,1,1),[1,5,10],RotorSpeed(1)/2/pi,'titletext','Orbit Nodes 1,5 and 10');
subplot(1,2,2)
plotOrbit(Mode_sort(:,1,15),[1,5,10],RotorSpeed(15)/2/pi,'titletext','Orbit Node  1,5 and 10');

%% response
RotorSpeed = (10:10:3000)*pi/30;
response = FRF(rotor,RotorSpeed);
plotFRF(RotorSpeed*30/pi,response,[10.1 21.1])

%% critical models comparison
num_crit_plot = 6; 

figure('Name', 'Critical Modes', 'NumberTitle', 'off')

for ii = 1:num_crit_plot
    subplot(2, 3, ii)
    freq_hz = critical_speeds(ii) / (2*pi);
    speed_rpm = critical_speeds(ii) * 30/pi;
    plotMode(rotor, mode_shapes(:, ii), ...
        'number', ii, ...
        'frequency', freq_hz, ...
        'rotorSpeed', speed_rpm);
    
    title(sprintf('Mode %d - %.2f rpm', ii, speed_rpm));
end
%% rappresentazione delle orbite
figure('Name','Orbits','NumberTitle', 'off')
subplot(1,2,1)
plotOrbit(mode_shapes(:,1),[1,5,10],critical_speeds(1)/2/pi,'titletext','Orbit Nodes 1,5 and 10');
subplot(1,2,2)
plotOrbit(mode_shapes(:,3),[1,5,10],critical_speeds(3)/2/pi,'titletext','Orbit Node  1,5 and 10');

%% root Locus
RotorSpeed = 10:10:3000;
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(1,:));

%% Stress Analysis
[nu_static,nu_fatigue,N_life] = stressAnalysis(rotor,response(:,100),350e6);
plotStressAnalysis(rotor,response(:,100),nu_static); 
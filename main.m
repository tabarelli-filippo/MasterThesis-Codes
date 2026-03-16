clear
close all
clc

addpath(genpath('auxFunc/'))

%% Generate the rotor
% shaft elements
coord_should = [0 20 34 42 52.9 93.9 108.3 160.3 176.5 184.7 185.7 187.7 193.9 233.3]'*1e-3; %[m] shaft
d_Ext1Shaft = [6 7 6.2 9.5 9.825 10 9.825 10 12.7 19.3 34.3 34.3 18.5]'* 2e-3; %[m] shaft shoulders
d_Ext2Shaft = [6 7 6.2 9.5 9.825 10 9.825 10 12.7 21.3 34.3 34.3 18.5]'* 2e-3; %[m] shaft shoulders
d_IntShaft = zeros(size(d_Ext1Shaft));

% shaft properties
rho = 7800; %[Kg/m^3]
E = 200e9; %[Pa]
poisson = 0.285; %[-]
G = E / 2 / (1 + poisson); %[Pa]

%% Disks
% compressor
z_compr = 50 *1e-3; %[m] 
m_compr = 2300*1e-3; %[Kg]
Id_compr = 3.5e6 * 1e-9; Ip_compr = 4e6 * 1e-9; %[Kg*m^2] 
d_ext_shaft_comp = 9.5*2e-3; %[m] % shaft diameter at compressor
d_int_shaft_comp = 0; %int shaft diameter at compressor

% turbine
z_turb = 214.71 *1e-3; %[m] 
m_turb = 3650*1e-3; %[Kg]
Id_turb = 4e6 * 1e-9; Ip_turb = 8e6 * 1e-9; %[Kg*m^2]
d_extShaft_turb = 18.5*2e-3; % ext shaft diameter at turbine
d_intShaft_turb = 0; %int shaft diameter at turbine

z_disks = [z_compr, z_turb]';

%% Bearings
% ball bearing
% n_ball = 16;
% d_ball = 3e-3
% fs = 150*9.8

d_m_bearing = 0.3e-3; % [m]
alpha = 20*pi/180; % [rad] 
delta = 0.5*d_m_bearing*tan(alpha);
z_bearings = [(93.9+108.3)/2-delta (160.3+176.5)/2+delta]'*1e-3; %[m] %effective bearing position
d_extShaft_bear = 10*2e-3; %[m] % ext shaft diameter at location
d_intShaft_bear = 0; %[m] % int shaft diameter at location

%% riodinamento coordinate e diametri
[z_coords, idx_sort] = sort([coord_should;z_disks;z_bearings]);

numEl = length(z_coords)-1;
zero_vec = zeros(numEl,1); %vettore ausiliario

d_extShaft_augm1 = [d_Ext1Shaft; d_ext_shaft_comp;d_extShaft_turb; d_extShaft_bear; d_extShaft_bear];
d_extShaft_augm2 = [d_Ext2Shaft; d_ext_shaft_comp;d_extShaft_turb; d_extShaft_bear; d_extShaft_bear];
d_intShaft_augm = [d_IntShaft; d_int_shaft_comp;d_intShaft_turb; d_intShaft_bear; d_intShaft_bear];

% definizione diametri ai nodi degli elementi albero
d_shaftEl = [zero_vec,d_extShaft_augm1,zero_vec,d_extShaft_augm2]; 
d_shaftEl = d_shaftEl(idx_sort(2:end)-1,:);
%% definizione elementi albero
% %numero di nodi intermedi aggiuntivi per elemento
n_intNodes = 2*ones(numEl,1);

% discretizzazione degli elementi infittita agli estremi dei nodi per elemento
isGraded = zero_vec; 
% tipo di elemento albero
shaftElType = 3*ones(numEl,1);
shaftElType(13) = 23; % elemento taper
% proprietà elemento albero
shaftElProperties = [rho*ones(numEl,1),E*ones(numEl,1),poisson*ones(numEl,1)];

%% Mesh Generator
rotor = meshGenerator(z_coords,d_shaftEl,n_intNodes,isGraded,shaftElType,shaftElProperties);

%% Disks and bearings definition
z_nodes_mesh = [rotor.nodes.node];
z_coord_mesh = [rotor.nodes.coord];

idx_compr = find(z_coord_mesh == z_compr,1);
idx_turb = find(z_coord_mesh == z_turb,1);
rotor.disk(1) = struct('type', 2, 'node', z_nodes_mesh(idx_compr), 'mass', m_compr, 'Id', Id_compr, 'Ip', Ip_compr);
rotor.disk(2) = struct('type', 2, 'node', z_nodes_mesh(idx_turb), 'mass', m_turb, 'Id', Id_turb, 'Ip', Ip_turb);

idx_bear1 = find(z_coord_mesh == z_bearings(1),1);
idx_bear2 = find(z_coord_mesh == z_bearings(2),1);
rotor.bearing(1) = struct('type', 1, 'node', z_nodes_mesh(idx_bear1));
rotor.bearing(2) = struct('type', 1, 'node', z_nodes_mesh(idx_bear2));

%% Draw the rotor
figureRotor(rotor);
mass = massRotor(rotor);
W = mass * 9.8061 / 2;

%% Campbell Diagram
rotorSpeed = (0:1000:10000)*pi/30; %[rad/s]

[natPuls,Mode,kappa] = charRoots(rotor,rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

NatPuls_sort0 = NatPuls_sort(1:1:7,:);
kappa_sort0 = kappa_sort(:,1:1:7,:);

NX = [1,2];
isDamped = 1;
plotCampbell(rotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);

%% plot Modes
figure('Name','Modes','NumberTitle', 'off')
subplot(1,3,1)
plotMode(rotor,Mode_sort(:,1,5),'number',1,'frequency',abs(imag(NatPuls_sort0(1,5)))/2/pi,'rotorSpeed',rotorSpeed(5)*30/pi);
subplot(1,3,2)
plotMode(rotor,Mode_sort(:,3,5),'number',3,'frequency',abs(imag(NatPuls_sort0(3,5)))/2/pi,'rotorSpeed',rotorSpeed(5)*30/pi);
subplot(1,3,3)
plotMode(rotor,Mode_sort(:,5,5),'number',5,'frequency',abs(imag(NatPuls_sort0(5,5)))/2/pi,'rotorSpeed',rotorSpeed(5)*30/pi);
%% plot Root Locus
plotRootLocus(rotorSpeed*30/pi,NatPuls_sort0(1,:));

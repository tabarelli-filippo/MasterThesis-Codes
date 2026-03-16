clear
close all

addpath(genpath('auxFunc/'))
%% Validating model doi:10.1016/j.measurement.2012.01.032
% Authors modelled shaft mass matrix as lumped, not as consistend, which
% leads to lower resonance frequencies. Correct formulation requires
% consistent matrix, as it is evaluated in this script. Natural frequencies
% here computed result in higher values. 

% Only forward modes are plotted in the paper, backward modes are drawn as 
% negative. No distinction in modes is made here, all critical speeds are 
% evaluated.

z_coords = [0 30 98.3 132.5 162.5 250.5 338.5 363.5 393.5 433.5 453.5]'*1e-3; %[m] shaft
d_ExtShaft = [30 120 120 63 68 68 150 150 40 30]'* 1e-3; %[m] shaft
d_IntShaft = [15 0 0 55 60 60 0 0 30 20]'* 1e-3;%[m] shaft

numEl = length(d_IntShaft);

% shaft properties
rho = [7810 0 0 7810 7810 7810 0 0 7810 7810]'; %[Kg/m^3]
%rho = 7810*ones(numEl,1);
E = 210e9; %[Pa]

%% Disks
% compressor
z_compr = 98.3*1e-3; %[m] 
m_compr = 20.81; %[Kg]
Id_compr = 0.174; Ip_compr = 0.285; %[Kg*m^2] 
d_ext_shaft_comp = 0; %[m] % shaft diameter at compressor
d_int_shaft_comp = 0; %int shaft diameter at compressor
epsilon = 1e-6; %[m]
% turbine
z_turb = 363.5e-3; %[m] 
m_turb = 18.2; %[Kg]
Id_turb = 0.142; Ip_turb = 0.269; %[Kg*m^2]
d_extShaft_turb = 0; % ext shaft diameter at turbine
d_intShaft_turb = 0; %int shaft diameter at turbine

z_disks = [z_compr, z_turb]';

%% Stiff Bearings
z_bearings = [0 453.5]'*1e-3; %[m]

%% Shaft elements definition
zero_vec = zeros(numEl,1);

% no internal discretization
n_intNodes = zero_vec;
isGraded = zero_vec;

% element type
shaftElType = 1*ones(numEl,1); %type 1 shaft element

% material shaft properties
shaftElProperties = [rho,E*ones(numEl,1),zero_vec];

% shaft diameters
d_shaftEl = [d_IntShaft, d_ExtShaft, d_IntShaft, d_ExtShaft];
%% Shaft Mesh
rotor = meshGenerator(z_coords,d_shaftEl,n_intNodes,isGraded,shaftElType,shaftElProperties);

%% Disks and bearings definition
z_nodes_mesh = [rotor.nodes.node];
z_coord_mesh = [rotor.nodes.coord];

idx_compr = find(z_coord_mesh == z_compr,1);
idx_turb = find(z_coord_mesh == z_turb,1);
rotor.disk(1) = struct('type', 2, 'node', z_nodes_mesh(idx_compr), 'mass', m_compr, 'Id', Id_compr, 'Ip', Ip_compr);
rotor.disk(2) = struct('type', 2, 'node', z_nodes_mesh(idx_turb), 'mass', m_turb, 'Id', Id_turb, 'Ip', Ip_turb);

rotor.disk(1).epsilon = epsilon;
rotor.disk(2).epsilon = 0;
rotor.forcing(1).type = 1;

idx_bear1 = find(z_coord_mesh == z_bearings(1),1);
idx_bear2 = find(z_coord_mesh == z_bearings(2),1);
rotor.bearing(1) = struct('type', 1, 'node', z_nodes_mesh(idx_bear1));
rotor.bearing(2) = struct('type', 1, 'node', z_nodes_mesh(idx_bear2));

%% Draw rotor
figureRotor(rotor);

%% Campbell Diagram
rotorSpeed = (0:1e3:1e5); %[rad/s]

[natPuls,Mode,kappa] = charRoots(rotor,rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

selectModes = 1:12;
NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0 = kappa_sort(:,selectModes,:);

NX = 1;
isDamped = 0;
plotCampbell(rotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);
%% Critical Speeds and critical modes
method = 2;
[criticalSpeeds,modeCrits] = critSpeeds(rotor,NX,isDamped,method,'num_crit',9);

% modes
figure('Name','Modes','NumberTitle', 'off')
selectModes = [2, 4, 8];
for jj = 1:length(selectModes)
    subplot(1,3,jj)
    kk = selectModes(jj);
    plotMode(rotor,modeCrits(:,kk),'number',kk,'frequency',abs(imag(NatPuls_sort0(1,kk)))/2/pi,'rotorSpeed',criticalSpeeds(kk)*30/pi);
    fprintf('Critical Speed %d : %.3f [rpm]\n',kk,criticalSpeeds(kk)*30/pi);
end

%% unbalancing
rotorSpeed = 0:10:6000;
response = FRF(rotor,rotorSpeed);
plotFRF(rotorSpeed*30/pi,response,[3.1 3.2])

%% plot displacement
% 2nd critical speed
idx_crit = find(rotorSpeed == 2310,1);
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

% 4th critical speed 
idx_crit = find(rotorSpeed == 3750,1);
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

% 6th critical speed 
idx_crit = find(rotorSpeed == 5640,1);
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

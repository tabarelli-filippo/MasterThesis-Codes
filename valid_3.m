clear
close all

addpath(genpath('auxFunc/'))
%% Validating model http://dx.doi.org/10.1016/j.measurement.2014.03.010
Element_length = [0.0071; 0.0071; 0.0053; 0.0054; 0.005; 0.0061; 0.0068; 0.0069; 0.0069;...
     0.0072; 0.010; 0.010; 0.010; 0.010; 0.0084; 0.0084; 0.0069; 0.0057;...
     0.008; 0.0101; 0.0113; 0.008];

numEl = length(Element_length);
coord_should = cumsum(Element_length);

% external diameter
d_ExtShaft = [0.007; 0.007; 0.053; 0.053; 0.012; 0.012; 0.0132; 0.0144; 0.0158; ...
      0.017; 0.018; 0.018; 0.018; 0.018; 0.0165; 0.015; 0.0136; 0.010; ... 
      0.0055; 0.060; 0.030; 0.0065];
% internal diameter
d_IntShaft = zeros(22, 1);

% material properties
rho = [2770; 2770; 0; 0; 2770; 2770; 2770; 2770; 2770; 2770; 2770; 2770;...
       2770; 2770; 2770; 2770; 2770; 2770; 2770; 0; 0; 2770];

E = [72; 72; 200; 200; 72; 72; 72; 72; 72; 72;...
     72; 72; 72; 72; 72; 72; 72; 72; 72; 96; 96; 72]*1e9;

poisson = [0.334; 0.334; 0.295; 0.295; 0.334; 0.334; 0.334; 0.334; 0.334; 0.334;...
    0.334; 0.334; 0.334; 0.334; 0.334; 0.334; 0.334; 0.334; 0.334; 0.34; 0.34; 0.334;];


%% Disks
% turbine
z_turb = coord_should(3); %[m] 
m_turb = 74.79e-3; %[Kg]
Id_turb = 1.32e-5; Ip_turb = 2.61e-5; %[Kg*m^2]
d_extShaft_turb = 53e-3; % ext shaft diameter at turbine
d_intShaft_turb = 0; %int shaft diameter at turbine
epsilon = 0.07479e-6 / m_turb;

% compressor
z_compr = coord_should(20); %[m] 
m_compr = 84.42e-3; %[Kg]
Id_compr = 1.28e-5; Ip_compr = 1.91e-5; %[Kg*m^2] 
d_ext_shaft_comp = 30e-3; %[m] % shaft diameter at compressor
d_int_shaft_comp = 0; %int shaft diameter at compressor

z_disks = [z_turb z_compr]';
%% bearings
z_bearings = [coord_should(5) coord_should(17)]; % [m]
d_extShaft_bear = [25.4 38.1]'* 1e-3;  % [m]

%% Shaft Elements
z_coords = [0; coord_should];
zero_vec = zeros(numEl,1);

% no internal discretization
n_intNodes = zero_vec;
isGraded = zero_vec;

% element type
shaftElType = 3*ones(numEl,1); %type 3 - Timoshenko Hutchinson

% material shaft properties
shaftElProperties = [rho,E,poisson];

% shaft diameters
d_shaftEl = [d_IntShaft, d_ExtShaft, d_IntShaft, d_ExtShaft];

rotor = meshGenerator(z_coords,d_shaftEl,n_intNodes,isGraded,shaftElType,shaftElProperties);

%% Disks & Bearings
z_nodes_mesh = [rotor.nodes.node];
z_coord_mesh = [rotor.nodes.coord];

[~, idx_turb_mesh] = find(z_coord_mesh == z_turb,1);  
[~, idx_compr_mesh] = find(z_coord_mesh == z_compr,1);
[~, idx_bear1_mesh] = find(z_coord_mesh == z_bearings(1),1);
[~, idx_bear2_mesh] = find(z_coord_mesh == z_bearings(2),1);

rotor.disk(1) = struct('type', 2, 'node', z_nodes_mesh(idx_turb_mesh), 'mass', m_turb, 'Id', Id_turb, 'Ip', Ip_turb);
rotor.disk(2) = struct('type', 2, 'node', z_nodes_mesh(idx_compr_mesh), 'mass', m_compr, 'Id', Id_compr, 'Ip', Ip_compr);

rotor.disk(1).epsilon = epsilon; 
rotor.disk(2).epsilon = 0;

rotor.forcing(1).type = 1;

% rigid bearings
rotor.bearing(1) = struct('type', 1, 'node', z_nodes_mesh(idx_bear1_mesh));
rotor.bearing(2) = struct('type', 1, 'node', z_nodes_mesh(idx_bear2_mesh));

% k = 5e8; %[N/m]
% c = 0; %[Ns/m]
% rotor.bearing(1) = struct('type', 3, 'node', z_nodes_mesh(idx_bear1_mesh),'kx',k,'ky',k,'cx',c,'cy',c);
% rotor.bearing(2) = struct('type', 3, 'node', z_nodes_mesh(idx_bear2_mesh),'kx',k,'ky',k,'cx',c,'cy',c);
%% rotor figure
figureRotor(rotor);
mass = massRotor(rotor);
%% Campbell Diagram
rotorSpeed = (0:1.5e4:1.5e5)*pi/30; %[rad/s]

[natPuls,Mode,kappa] = charRoots(rotor,rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

selectModes = 1:6;
%selectModes = [1,3,6];

NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0 = kappa_sort(:,selectModes,:);

NX = 1;
isDamped = 0;
plotCampbell(rotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);

%% Critical Speeds
method = 1;
[criticalSpeeds,modeCrits] = critSpeeds(rotor,NX,isDamped,method);

for jj = 1:length(selectModes)
    fprintf('Critical Speed %d : %.3f [rpm]\n',jj,criticalSpeeds(jj)*30/pi);
end

% modes
figure('Name','Modes','NumberTitle', 'off')
selectModes = [1, 2, 3];
for jj = 1:length(selectModes)
    subplot(1,3,jj)
    kk = selectModes(jj);
    plotMode(rotor,modeCrits(:,kk),'number',kk,'frequency',abs(imag(NatPuls_sort0(1,kk)))/2/pi,'rotorSpeed',criticalSpeeds(kk)*30/pi);
    fprintf('Critical Speed %d : %.3f [rpm]\n',kk,criticalSpeeds(kk)*30/pi);
end

%% Unbalancing response
rotorSpeed = (0:1e2:1.5e5)*pi/30; 
response = FRF(rotor,rotorSpeed);
plotFRF(rotorSpeed*30/pi,response,[4.1 4.2 5.1 7.2 8.1])

%% plot displacement

% 2nd critical speed

[~, idx_crit] = min(abs(rotorSpeed - 25000*pi/30));
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

% 4th critical speed 
[~, idx_crit] = min(abs(rotorSpeed - 45700*pi/30));
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)

% 6th critical speed 
[~, idx_crit] = min(abs(rotorSpeed - 79800*pi/30));
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)



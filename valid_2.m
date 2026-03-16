clear
close all

addpath(genpath('auxFunc/'))
%% Validating model doi.org/10.1016/j.measurement.2018.08.044
coord_should = [0 101.6 254 508 558.8 660.4 711.2 812.8]'*1e-3; %[m] shaft
d_ExtShaft = [25.4 50.8 76.2 50.8 38.1 44.4 50.8]'* 1e-3; %[m] shaft
d_IntShaft = 0;%[m] shaft

% shaft properties
rho = 7833.41; %[Kg/m^3]
E = 206.84e9; %[Pa]

%% Disks
z_disks = [812.8/2 (812.8/2+50.8) (812.8-50.8/2)]'*1e-3; %[m]
D_disk = [558.8 297.4 254]'*1e-3; % [m]
thick_disk = [50.8 50.8 50.8]'*1e-3; % [m]
d_ext_shaft_disks = [76.2 76.2 50.8]'* 1e-3; %[m]

%% Bearings
z_bearings = [81.28 640.08]'* 1e-3; % [m]
d_extShaft_bear = [25.4 38.1]'* 1e-3;  % [m]

%% Shaft elements definition
[z_coords, idx_sort] = sort([coord_should;z_disks;z_bearings]);

numEl = length(z_coords)-1;
zero_vec = zeros(numEl,1); %vettore ausiliario

d_extShaft_augm = [d_ExtShaft; d_ext_shaft_disks; d_extShaft_bear;];

% definizione diametri ai nodi degli elementi albero
d_shaftEl = [zero_vec,d_extShaft_augm,zero_vec,d_extShaft_augm]; 
d_shaftEl = d_shaftEl(idx_sort(2:end)-1,:);

% no internal discretization
n_intNodes = zero_vec;
n_intNodes([9,12]) = 2;
n_intNodes([5,6,7,10,11]) = 3;
n_intNodes([1,3,4,8]) = 5;
isGraded = zero_vec;

% element type
shaftElType = 1*ones(numEl,1); %type 1 shaft element

% material shaft properties
shaftElProperties = [rho*ones(numEl,1),E*ones(numEl,1),zero_vec];

%% Shaft Mesh
rotor = meshGenerator(z_coords,d_shaftEl,n_intNodes,isGraded,shaftElType,shaftElProperties);

%% Disks
rotor.disk(3) = struct('type', [], 'node', [], 'thick', [], 'D_int', ...
    [], 'D_ext', [], 'rho',[],'m0',[],'delta',[]);

[rotor.disk.type] = deal(1);
[rotor.disk.rho] = deal(rho);

d_ext_shaft_disks_cell = num2cell(d_ext_shaft_disks);
[rotor.disk.D_int] = d_ext_shaft_disks_cell{:}; %shaft diameters
D_disk_cell = num2cell(D_disk);
[rotor.disk.D_ext] = D_disk_cell{:}; % disks diameters
thick_disk_cell = num2cell(thick_disk);
[rotor.disk.thick] = thick_disk_cell{:}; %thick diameters

%disk nodes
z_nodes_mesh = [rotor.nodes.node];
z_coord_mesh = [rotor.nodes.coord];

idx_d1 = find(z_coord_mesh == z_disks(1),1);
idx_d2 = find(z_coord_mesh == z_disks(2),1);
idx_d3 = find(z_coord_mesh == z_disks(3),1);

node_disks = [idx_d1, idx_d2, idx_d3];
node_disks_cell = num2cell(node_disks);
[rotor.disk.node] = node_disks_cell{:};

%unbalancing
m0_disks = [0 2.54e-4 5.08e-5]./D_disk'*2; %[Kg]
m0_phase = [0 pi/2 0]; %[rad]
m0_disks_cell = num2cell(m0_disks);
[rotor.disk.m0] = m0_disks_cell{:};
m0_phase_cell = num2cell(m0_phase);
[rotor.disk.delta] = m0_phase_cell{:};

%% Bearings
idx_b1 = find(z_coord_mesh == z_bearings(1),1);
idx_b2 = find(z_coord_mesh == z_bearings(2),1);

k = 5.534e6; %[N/m]
c = 1.7513e3; %[Ns/m]

rotor.bearing(1) = struct('type', 3, 'node', z_nodes_mesh(idx_b1),'kx',k,'ky',k,'cx',c,'cy',c);
rotor.bearing(2) = struct('type', 3, 'node', z_nodes_mesh(idx_b2),'kx',k,'ky',k,'cx',c,'cy',c);
%% forcing
rotor.forcing(1).type = 1;
%% figure 
figureRotor(rotor);

%% Campbell Diagram
rotorSpeed = (0:1000:6000)*pi/30; %[rad/s]

[natPuls,Mode,kappa] = charRoots(rotor,rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

selectModes = 1:1:6;
NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0 = kappa_sort(:,selectModes,:);

NX = 1;
isDamped = 1;
plotCampbell(rotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);

%% Modes
% modes
figure('Name','Modes','NumberTitle', 'off')

for jj = 1:length(selectModes)
    subplot(2,3,jj)
    kk = selectModes(jj);
    plotMode(rotor,Mode_sort(:,kk,2),'number',kk,'frequency',abs(NatPuls_sort0(1,kk))/2/pi,'rotorSpeed',rotorSpeed(2)*30/pi);
end

%% unbalancing
rotorSpeed = (0:10:6000)*pi/30;
response = FRF(rotor,rotorSpeed);

plotFRF(rotorSpeed*30/pi,response,[24.1 24.2])

%% Critical Speeds vs Stiffness
stiff = logspace(log10(3e6), log10(2e9), 7);%[N/m]
nn = length(stiff);
NX = 1;
criticalSpeeds = zeros(2,nn);

for jj = 1:nn
    k = stiff(jj);
    rotor.bearing(1) = struct('type', 3, 'node', z_nodes_mesh(idx_b1),'kx',k,'ky',k,'cx',c,'cy',c);
    rotor.bearing(2) = struct('type', 3, 'node', z_nodes_mesh(idx_b2),'kx',k,'ky',k,'cx',c,'cy',c);
    
    method = 2;
    [Speeds,~] = critSpeeds(rotor,NX,isDamped,method,'num_crit',4);
    criticalSpeeds(:,jj) = Speeds([2,4]);
end

figure('Name','Critical Speeds Map','NumberTitle', 'off', 'Color', 'w');
clf;

hold on;

loglog(stiff, criticalSpeeds(1,:)*30/pi, 'b-o', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', '1st Mode');
loglog(stiff, criticalSpeeds(2,:)*30/pi, 'r-s', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', '2nd Mode');
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log');
grid on; grid minor;
xlabel('Bearing Stiffness k [N/m]');
ylabel('Critical Speed [RPM]'); 
title('Critical Speed Map vs Bearing Stiffness');
legend('show', 'Location', 'best');
xlim([min(stiff)*0.8, max(stiff)*1.2]);
ylim([1e3,1.5e4])

%% Stability Analysis
rotor = rmfield(rotor, 'bearing');
% bearing characteristics
rotorSpeed_points = [1 1000 3000 5000 7000 10000]*pi/30; %[rad/s]

zero_points = zeros(1,6);
kxx_points = [1.58e6 1.40e6 1.23e6 1.05e6 8.76e5 7.00e5]; % [N/m] stiffness coeff points
kyy_points = [1.58e6 1.40e6 1.23e6 1.05e6 8.76e5 7.00e5]; % [N/m] stiffness coeff points
cxx_points = [1.75e3 5.25e3 8.76e3 1.23e4 1.58e4 1.93e4]; % [Ns/m] damping coeff points
cyy_points = [1.75e3 5.25e3 8.76e3 1.23e4 1.58e4 1.93e4]; % [Ns/m] damping coeff points

rotor.bearing(1) = struct('type', 13, 'node', z_nodes_mesh(idx_b1),...
    'Speed_points',rotorSpeed_points,...
    'kxx_points',kxx_points,'kxy_points',zero_points,...
    'kyx_points',zero_points,'kyy_points',kyy_points, ...
    'cxx_points',cxx_points,'cxy_points',zero_points,...
    'cyx_points',zero_points,'cyy_points',cyy_points);


rotor.bearing(2) = struct('type', 13, 'node', z_nodes_mesh(idx_b2),...
    'Speed_points',rotorSpeed_points,...
    'kxx_points',kxx_points,'kxy_points',zero_points,...
    'kyx_points',zero_points,'kyy_points',kyy_points, ...
    'cxx_points',cxx_points,'cxy_points',zero_points,...
    'cyx_points',zero_points,'cyy_points',cyy_points);

% natural frequencies evaluation
RotorSpeed = (10:100:5000)*pi/30;
[natPuls,Mode,kappa] = charRoots(rotor,RotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

selectModes = 1:1:6;
NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0 = kappa_sort(:,selectModes,:);

%plotting root locus
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(1,:));
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(2,:));
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(3,:));
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(4,:));
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(5,:));
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(6,:));

NX = 1;
isDamped = 1;
plotCampbell(RotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);
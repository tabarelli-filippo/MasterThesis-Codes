clear
close all

addpath(genpath('auxFunc/'))
%% Validating model 10.2514/1.B34104

z_coords = [0; 0.04; 0.075; 0.087; 0.1; 0.157; 0.21; 0.253; 0.292; 0.36;...
    0.63; 0.642; 0.7; 0.735; 0.765; 0.78]; %[m]
numEl = length(z_coords) - 1;

d_ExtShaft = [0.015; 0.015; 0.040; 0.040; 0.040; 0.040; 0.040; ...
    0.040; 0.040; 0.034; 0.034; 0.015; 0.015; 0.0275; 0.0275]*2; %[m];

d_IntShaft = [0; 0; 0; 0.029; 0.029; 0.029; 0.029; 0.029; 0.029;...
    0.029; 0; 0; 0; 0; 0]*2; %[m]

rho = 7850*ones(numEl,1); %[Kg/m^3]
E = 2.05e11*ones(numEl,1); %[Pa]
poisson = 0.29*ones(numEl,1); %[-]
beta = zeros(numEl,1);
%% Disks and bearings
z_disks = [z_coords(5:9) ; z_coords(15)]; %[m] 
rho_disks = [2830; 2830; 2870; 4430; 4430; 8000];  %[Kg/m^3]
thick_disks = [0.02; 0.02; 0.02; 0.02; 0.02; 0.03]; %[m]
Dext_disks = [0.19; 0.21; 0.23; 0.24; 0.25; 0.25]; %[m]
d_extShaft_disks = [0.08; 0.08; 0.08; 0.08; 0.08; 0.055;]; % ext shaft diameter at disk

m0 = [0; 0; 0; 0; 0; 1e-4/0.125]; %[Kg]

% bearings
z_bearings = [z_coords(2); z_coords(13)];

%% Shaft Elements
zero_vec = zeros(numEl,1);

% no internal discretization
n_intNodes = 0*ones(numEl,1);
n_intNodes(10) = 6;
isGraded = zeros(numEl,1);

% element type
shaftElType = 3*ones(numEl,1); %type 3 - Timoshenko Hutchinson

% material shaft properties
shaftElProperties = [rho, E, poisson];

% shaft diameters
d_shaftEl = [d_IntShaft, d_ExtShaft, d_IntShaft, d_ExtShaft];

% shaft load
shaftLoad = zeros(numEl,2);

rotor = meshGenerator(z_coords,d_shaftEl,n_intNodes,isGraded,shaftElType,shaftElProperties,shaftLoad,beta);

%% Disks
z_nodes_mesh = [rotor.nodes.node];
z_coord_mesh = [rotor.nodes.coord];

[found_flag, idx_disks_mesh] = ismember(z_disks, z_coord_mesh);
idx_disks_mesh = idx_disks_mesh(found_flag);

[~, idx_bear1_mesh] = find(z_coord_mesh == z_bearings(1),1);
[~, idx_bear2_mesh] = find(z_coord_mesh == z_bearings(2),1);

rotor.disk(6) = struct('type', [], 'node', [], 'thick', [], 'D_ext', [], ...
    'D_int', [], 'rho', [], 'm0',[]);

% disks type
[rotor.disk.type] = deal(1);
% disks node
node_disks_cell = num2cell(idx_disks_mesh);
[rotor.disk.node] = node_disks_cell{:};
% disks thick & diameters
thick_disks_cell = num2cell(thick_disks);
[rotor.disk.thick] = thick_disks_cell{:};

Dext_disks_cell = num2cell(Dext_disks);
[rotor.disk.D_ext] = Dext_disks_cell{:};

Dint_disks_cell = num2cell(d_extShaft_disks);
[rotor.disk.D_int] = Dint_disks_cell{:};
% density
rho_disks_cell = num2cell(rho_disks);
[rotor.disk.rho] = rho_disks_cell{:};
% unbalancing
m0_cell = num2cell(m0);
[rotor.disk.m0] = m0_cell{:};


%% bearings
rotorSpeed_points = linspace(0,500,7)*2*pi; % speed sampling points
% front bearing 
kxx1_points = 1.37e7*ones(1,7); %[N/m] stiffness coeff points
kxy1_points = zeros(1,7); % [N/m]stiffness coeff points
cxx1_points = [39.3 33.4 19.4 13.5 9.3 7.0 5.1]; % [Ns/m] stiffness coeff points
cxy1_points = zeros(1,7); % [Ns/m] stiffness coeff points

rotor.bearing(1) = struct('type', 13, 'node', z_nodes_mesh(idx_bear1_mesh),...
    'Speed_points',rotorSpeed_points,...
    'kxx_points',kxx1_points,'kxy_points',kxy1_points,...
    'kyx_points',kxy1_points,'kyy_points',kxx1_points, ...
    'cxx_points',cxx1_points,'cxy_points',cxy1_points,...
    'cyx_points',cxy1_points,'cyy_points',cxx1_points);

% rear bearing 
kxx2_points = [2.94 2.94 5.30 7.16 8.44 10.5 12.5]*1e6; %[N/m] stiffness coeff points
kxy2_points = zeros(1,7); % [N/m]stiffness coeff points
cxx2_points = [9719.7 9719.7 7855.7 6055.9 4216.1 3352.2 2441.2]; % [Ns/m] stiffness coeff points
cxy2_points = zeros(1,7); % [Ns/m] stiffness coeff points

rotor.bearing(2) = struct('type', 13, 'node', z_nodes_mesh(idx_bear2_mesh),...
    'Speed_points',rotorSpeed_points,...
    'kxx_points',kxx2_points,'kxy_points',kxy2_points,...
    'kyx_points',kxy2_points,'kyy_points',kxx2_points, ...
    'cxx_points',cxx2_points,'cxy_points',cxy2_points,...
    'cyx_points',cxy2_points,'cyy_points',cxx2_points);

rotor.forcing(1) = struct('type', 1);
%% rotor figure
figureRotor(rotor);

%% Campbell Diagram
rotorSpeed = linspace(0,500,10)*2*pi; %[rad/s]

[natPuls,Mode,kappa] = charRoots(rotor,rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

selectModes = 1:6;

NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0 = kappa_sort(:,selectModes,:);

NX = [1 0.5];
isDamped = 1;
plotCampbell(rotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);
%% Critical Modes
method = 2;
NX = 1;
[criticalSpeeds,modeCrits] = critSpeeds(rotor,NX,isDamped,method,"num_crit",6);

% modes
figure('Name','Modes','NumberTitle', 'off')

selectModes = 1:6;
for jj = 1:length(selectModes)
    subplot(2,3,jj)
    kk = selectModes(jj);
    plotMode(rotor,modeCrits(:,kk),'number',kk,'frequency',criticalSpeeds(kk)/2/pi,'rotorSpeed',criticalSpeeds(kk)*30/pi);
    fprintf('Critical Speed %d : %.3f [rpm]\n',kk,criticalSpeeds(kk)*30/pi);
end

%% Stability plot
% 1FW Mode Stabilty
plotRootLocus(rotorSpeed*30/pi,NatPuls_sort(2,:))
% 1FW Mode Stabilty
plotRootLocus(rotorSpeed*30/pi,NatPuls_sort(4,:))
% LogDecrement 
deltaLog = LogDecrement(NatPuls_sort);
plotLogDecrement(rotorSpeed/2/pi, NatPuls_sort(1,:))


%% Response to unbalancing
rotorSpeed = (0:0.1:500)*2*pi; 
response = FRF(rotor,rotorSpeed);
plotFRF(rotorSpeed*30/pi,response,[5.1,5.2])
plotFRF(rotorSpeed*30/pi,response,[21.1,21.2])

[~, idx_crit] = min(abs(rotorSpeed - criticalSpeeds(4)));
plotDisplacement(rotor, response(:,idx_crit), 'rotorSpeed', rotorSpeed(idx_crit)*30/pi)
%% RUNUP
alpha = [0.5*100*pi 0 0];
t_span = [0,10];
t_end = t_span(end);

numEl = length(z_nodes_mesh);

initPos = zeros(4*(numEl),1);
initVel = zeros(4*(numEl),1);

[time,q,qdot,speed] = runUp(rotor,t_span,alpha,initPos,initVel,nr = 12);

%% 2. Preparazione Dati
% Definisco i gradi di libertà per chiarezza
node_target = 21;
idx_u = (node_target-1)*4 + 1; % Spostamento orizzontale (x)
idx_v = (node_target-1)*4 + 2; % Spostamento verticale (y)

disp_u = q(:, idx_u);
disp_v = q(:, idx_v);
disp_total = sqrt(disp_u.^2 + disp_v.^2);
rpm_axis = speed * 30 / pi;

%% 3. Figura 1: Risposta nel Tempo (Transient Response)
fig1 = figure('Name','Runup_Time','NumberTitle', 'off', 'Color', 'w');
plot(time, disp_v, 'Color', [0 0.4470 0.7410], 'LineWidth',1.5) % Blu standard MATLAB
grid on; grid minor;
box on; % Chiude il riquadro del grafico

% Etichette e Titoli con unità di misura corrette
xlabel('Time t [s]');
ylabel('Vertical Displacement v [m]');
title(['Transient Response at Node ' num2str(node_target)]);

% Limiti (opzionale: adatta se necessario)
xlim([0, max(time)]);

%% 4. Figura 2: Risposta in Frequenza (Unbalance Response vs RPM)
fig2 = figure('Name','Runup_RPM','NumberTitle', 'off', 'Color', 'w');
plot(rpm_axis, disp_total, 'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5) % Rosso arancio
grid on; grid minor;
box on;

xlabel('Rotor Speed \Omega [rpm]');
ylabel('Total Displacement Amplitude r [m]');
title(['Unbalance Response (Run-up) - Node ' num2str(node_target)]);

% Opzionale: Aggiungi linea per velocità critiche se le conosci
% xline(critical_speed_rpm, '--k', 'Critical Speed');


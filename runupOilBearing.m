clear
close all

addpath(genpath('auxFunc/'))
% recreating example based on doi:10.1016/j.apm.2004.09.003
%% Jeffcott rotor
de = 0.012; %[m]
Le = 0.145; %[m]
L = 0.010; %[m]
coord_should = [L; 4*Le+L; 4*Le+2*L]; 

rho = 7800; %[Kg/m^3]
E = 206.7e9; %[Pa]
beta = 25e-5;
poisson = 0.285; %[-]
G = E / 2 / (1 + poisson); %[Pa]

d_shaftEl = repmat([0 de 0 de],3,1);
z_coords = [0 ;coord_should];
n_intNodes = [1;3;1];
isGraded = zeros(3,1);
shaftElType = 3*ones(3,1);
shaftElProperties = repmat([rho E poisson],3,1);
shaftLoad = zeros(3,2);
intDamping = beta * ones(3,1);

rotor = meshGenerator(z_coords,d_shaftEl,n_intNodes,isGraded,shaftElType,shaftElProperties,shaftLoad,intDamping);
%% Disks
z_disk = L + 2*Le; %[m] 
Dext_disk = 0.095; %[m]
thick_disk = 0.043; %[m]
m0 = 2e-5/Dext_disk; %[Kg]

z_coord_mesh = [rotor.nodes.coord];
z_nodes_mesh = [rotor.nodes.node];
idx_disk = find(z_coord_mesh == z_disk,1);

rotor.disk(1) = struct('type', 1, 'node', z_nodes_mesh(idx_disk), 'D_ext', Dext_disk, 'D_int',...
    de, 'thick', thick_disk, 'rho', rho , 'm0', m0);

%% Bearings
rotor.bearing(2) = struct('type', [], 'node', []);

idx_b1 = find(z_coord_mesh == L/2,1);
idx_b2 = find(z_coord_mesh == 1.5*L + 4*Le ,1);
node_bears_cell = num2cell([z_nodes_mesh(idx_b1),z_nodes_mesh(idx_b2)]);

D = 0.031; %[m]
c = 90e-6; %[m]
eta = 0.04; %[Pa s] viscosità olio

[rotor.bearing.type] = deal(7.1);
[rotor.bearing.node] = node_bears_cell{:};
[rotor.bearing.L] = deal(L);
[rotor.bearing.D] = deal(D);
[rotor.bearing.c] = deal(c);
[rotor.bearing.eta] = deal(eta);

%% forcing
rotor.forcing(1).type = 1;
rotor.forcing(2).type = 4;

%% figure
figureRotor(rotor);

%% solving
alpha = [4 0.01 0];
t_span = [0,15];
t_end = t_span(end);

numEl = length(z_nodes_mesh);

initPos = zeros(4*(numEl),1);
initPos(1:4:end) = 1e-7;
initVel = zeros(4*(numEl),1);
initVel(1:4:end) = 1e-7;

[time,q,qdot,speed] = runUp(rotor,t_span,alpha,initPos,initVel,nr = 4);

% figure('Name','Runup','NumberTitle', 'off')
% plot(speed*30/pi,q(:,1))
% xlabel('Rotor Speed [rpm]')
% ylabel('u disp[m]')
% title('u displacement - node 1')
%% Post-Processing and Visualization

% 1. Identify Disk Node (Critical DoF)
% We plot the disk node as it experiences the maximum amplitude during resonance.
node_disk = rotor.disk(1).node; 
idx_disk_plot = find([rotor.nodes.node] == node_disk, 1);

% Assuming 4 DoFs per node: [x, y, theta_x, theta_y]
dof_x = 4*(idx_disk_plot-1) + 1;
dof_y = 4*(idx_disk_plot-1) + 2;

% Extract displacements
disp_x = q(:, dof_x);
disp_y = q(:, dof_y);

% Calculate Orbit Amplitude (Instantaneous Envelope)
% r(t) = sqrt(x(t)^2 + y(t)^2)
amplitude = sqrt(disp_x.^2 + disp_y.^2);

% Convert speed to RPM
rpm = speed * 30 / pi;

%% Figure 1: Run-up Response (Bode Plot - Magnitude)
figure('Name', 'Run-up Response', 'NumberTitle', 'off', 'Color', 'w');

plot(rpm, disp_x, 'DisplayName', 'Disp X', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.0);
hold on;
plot(rpm, disp_y, 'DisplayName', 'Disp Y', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.0);
plot(rpm, amplitude, 'k-', 'DisplayName', 'Total Amplitude (r)', 'LineWidth', 1.5);
hold off;

grid on; grid minor;
xlabel('Rotor Speed [RPM]', 'FontWeight', 'bold');
ylabel('Displacement [m]', 'FontWeight', 'bold');
title(['Run-up Response at Disk Node (Node ' num2str(idx_disk_plot) ')']);
legend('Location', 'best');
xlim([0, max(rpm)]);

%% Figure 2: Rotor Orbit (Whirl)
figure('Name', 'Orbit Analysis', 'NumberTitle', 'off', 'Color', 'w');

plot(disp_x, disp_y, 'k', 'LineWidth', 0.5);
hold on;
% Mark Start and End points for direction clarity
plot(disp_x(1), disp_y(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(disp_x(end), disp_y(end), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', 'End');
hold off;

grid on; grid minor;
axis square; 
xlabel('Displacement X [m]', 'FontWeight', 'bold');
ylabel('Displacement Y [m]', 'FontWeight', 'bold');
title(['Rotor Orbit (Whirl) at Disk Node (Node ' num2str(idx_disk_plot) ')']);
legend('Orbit Path', 'Start', 'End', 'Location', 'best');
% VALID_3  Model validation against published micro-turbocharger data.
%   Reference: http://dx.doi.org/10.1016/j.measurement.2014.03.010
%
%   Validates the FE model against a high-speed micro-turbocharger benchmark
%   with mixed aluminum/steel shaft and two rigid bearings. The model uses
%   Timoshenko beam elements with Hutchinson shear coefficients (type 3)
%   and element-by-element material properties.
%
%   ROTOR DESCRIPTION
%     Shaft : 22 elements of varying length and diameter (aluminum and
%             steel sections); some elements have rho = 0 (structural
%             filler) with appropriate non-zero E
%             Total length ≈ 186 mm
%     Disks : 2 rigid disks (type 2):
%               Turbine:    m = 74.79 g at z ≈ element 3 end
%               Compressor: m = 84.42 g at z ≈ element 20 end
%               Turbine unbalance: epsilon = 0.07479e-6 / m_turb
%     Bearings: 2 pinned (type 1) at elements 5 and 17 ends
%     Elements: Timoshenko type 3, no intermediate refinement
%
%   ANALYSIS PERFORMED
%     1. Mesh, disk and bearing assembly; rotor schematic, mass
%     2. Campbell diagram for 6 sorted modes, 1X line (undamped)
%     3. Critical speeds (direct method 1, 6 speeds) with mode shapes
%        for modes 1, 2, 3 (plotMode)
%     4. Unbalance response at 5 DOFs (FRF, plotFRF)
%     5. Deflected shape at 3 critical speed regions (plotDisplacement)
%
% SEE ALSO
%   valid_1, valid_2, valid_4, charRoots, critSpeeds, FRF

clear
close all
addpath(genpath('auxFunc/'))

%% Shaft element lengths and diameters
Element_length = [0.0071;0.0071;0.0053;0.0054;0.005;0.0061;0.0068;0.0069;0.0069;...
    0.0072;0.010;0.010;0.010;0.010;0.0084;0.0084;0.0069;0.0057;...
    0.008;0.0101;0.0113;0.008];
numEl      = length(Element_length);
coord_should = cumsum(Element_length);

d_ExtShaft = [0.007;0.007;0.053;0.053;0.012;0.012;0.0132;0.0144;0.0158;...
    0.017;0.018;0.018;0.018;0.018;0.0165;0.015;0.0136;0.010;...
    0.0055;0.060;0.030;0.0065];
d_IntShaft = zeros(22,1);

rho     = [2770;2770;0;0;2770;2770;2770;2770;2770;2770;2770;2770;...
           2770;2770;2770;2770;2770;2770;2770;0;0;2770];
E       = [72;72;200;200;72;72;72;72;72;72;72;72;...
           72;72;72;72;72;72;72;96;96;72]*1e9;
poisson = [0.334;0.334;0.295;0.295;0.334;0.334;0.334;0.334;0.334;0.334;...
    0.334;0.334;0.334;0.334;0.334;0.334;0.334;0.334;0.334;0.34;0.34;0.334];

%% Disk and bearing parameters
z_turb   = coord_should(3);  m_turb  = 74.79e-3; Id_turb  = 1.32e-5; Ip_turb  = 2.61e-5;
z_compr  = coord_should(20); m_compr = 84.42e-3; Id_compr = 1.28e-5; Ip_compr = 1.91e-5;
epsilon  = 0.07479e-6 / m_turb;
z_bearings = [coord_should(5), coord_should(17)];

%% Mesh (no refinement, Timoshenko type 3)
z_coords = [0; coord_should];
zero_vec = zeros(numEl,1);
d_shaftEl = [d_IntShaft, d_ExtShaft, d_IntShaft, d_ExtShaft];
shaftElProperties = [rho, E, poisson];

rotor = meshGenerator(z_coords, d_shaftEl, zero_vec, zero_vec, 3*ones(numEl,1), shaftElProperties);

%% Disk and bearing assembly
z_nodes_mesh = [rotor.nodes.node]; z_coord_mesh = [rotor.nodes.coord];
idx_turb  = find(z_coord_mesh==z_turb, 1);
idx_compr = find(z_coord_mesh==z_compr, 1);
idx_b1 = find(z_coord_mesh==z_bearings(1), 1);
idx_b2 = find(z_coord_mesh==z_bearings(2), 1);

rotor.disk(1) = struct('type',2,'node',z_nodes_mesh(idx_turb), 'mass',m_turb, 'Id',Id_turb, 'Ip',Ip_turb);
rotor.disk(2) = struct('type',2,'node',z_nodes_mesh(idx_compr),'mass',m_compr,'Id',Id_compr,'Ip',Ip_compr);
rotor.disk(1).epsilon = epsilon; rotor.disk(2).epsilon = 0;
rotor.forcing(1).type = 1;
rotor.bearing(1) = struct('type',1,'node',z_nodes_mesh(idx_b1));
rotor.bearing(2) = struct('type',1,'node',z_nodes_mesh(idx_b2));

figureRotor(rotor);
mass = massRotor(rotor);
fprintf('Total rotor mass: %.4f kg\n', mass);

%% Campbell diagram
rotorSpeed = (0:1.5e4:1.5e5)*pi/30;
[natPuls, Mode, kappa] = charRoots(rotor, rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);
selectModes   = 1:6;
NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0   = kappa_sort(:,selectModes,:);
NX = 1; isDamped = false;
plotCampbell(rotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort0);

%% Critical speeds (direct method) and mode shapes
method = 1;
[criticalSpeeds, modeCrits] = critSpeeds(rotor, NX, isDamped, method);
for jj = 1:length(selectModes)
    fprintf('Critical Speed %d: %.3f rpm\n', jj, criticalSpeeds(jj)*30/pi);
end

figure('Name','Modes','NumberTitle','off')
selectModes_plot = [1, 2, 3];
for jj = 1:length(selectModes_plot)
    subplot(1,3,jj); kk = selectModes_plot(jj);
    plotMode(rotor, modeCrits(:,kk), 'number', kk, ...
        'frequency', abs(imag(NatPuls_sort0(1,kk)))/2/pi, ...
        'rotorSpeed', criticalSpeeds(kk)*30/pi);
end

%% Unbalance response
rotorSpeed = (0:1e2:1.5e5)*pi/30;
response   = FRF(rotor, rotorSpeed);
plotFRF(rotorSpeed*30/pi, response, [4.1 4.2 5.1 7.2 8.1])

%% Deflected shapes at selected critical speeds
[~,idx] = min(abs(rotorSpeed - 25000*pi/30));
plotDisplacement(rotor, response(:,idx), 'rotorSpeed', rotorSpeed(idx)*30/pi)

[~,idx] = min(abs(rotorSpeed - 45700*pi/30));
plotDisplacement(rotor, response(:,idx), 'rotorSpeed', rotorSpeed(idx)*30/pi)

[~,idx] = min(abs(rotorSpeed - 79800*pi/30));
plotDisplacement(rotor, response(:,idx), 'rotorSpeed', rotorSpeed(idx)*30/pi)

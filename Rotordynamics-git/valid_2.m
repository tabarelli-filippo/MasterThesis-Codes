% VALID_2  Model validation against published rotor dynamics data.
%   Reference: doi.org/10.1016/j.measurement.2018.08.044
%
%   Validates the FE model against the three-disk rotor benchmark from the
%   reference paper. Includes a full stability analysis using speed-dependent
%   bearing properties (type 13) after the frequency-domain validation.
%   Also generates a critical speed map as a function of isotropic bearing
%   stiffness.
%
%   ROTOR DESCRIPTION
%     Shaft  : 7 coarse sections, variable diameter (25.4–76.2 mm, solid)
%              refined mesh with 3–5 intermediate nodes per element
%     Disks  : 3 disks at z = L/2, L/2+50.8 mm, L-25.4 mm
%              D_ext = 558.8, 297.4, 254 mm; thick = 50.8 mm each
%              Unbalances: disk 2 (m0 = 2.54e-4 kg, phi = 90°),
%                          disk 3 (m0 = 5.08e-5 kg)
%     Bearings: 2 isotropic elastic+damping (type 3) initially at
%               z = 81.28 mm and z = 640.08 mm; k = 5.534e6 N/m,
%               c = 1751.3 N·s/m
%
%   ANALYSIS PERFORMED (Part 1 — Frequency Domain)
%     1. Mesh generation and figure
%     2. Campbell diagram (sorted, 6 modes, 1X, isDamped=true)
%     3. Mode shapes at speed index 2
%     4. Unbalance response at nodes 24.1 and 24.2 (FRF, plotFRF)
%     5. Critical speed map vs bearing stiffness (7 stiffness values,
%        2 modes tracked) — log-log plot
%
%   ANALYSIS PERFORMED (Part 2 — Stability)
%     Bearing replaced with speed-dependent type 13 (tabulated K and C).
%     6. Root locus for 6 modes (plotRootLocus)
%     7. Campbell diagram with speed-dependent bearings
%
% SEE ALSO
%   valid_1, valid_3, valid_4, charRoots, critSpeeds, plotRootLocus

clear
close all
addpath(genpath('auxFunc/'))

%% Shaft geometry
coord_should = [0 101.6 254 508 558.8 660.4 711.2 812.8]'*1e-3;
d_ExtShaft   = [25.4 50.8 76.2 50.8 38.1 44.4 50.8]'*1e-3;
rho = 7833.41; E = 206.84e9;

%% Disk parameters
z_disks           = [812.8/2 (812.8/2+50.8) (812.8-50.8/2)]'*1e-3;
D_disk            = [558.8 297.4 254]'*1e-3;
thick_disk        = [50.8 50.8 50.8]'*1e-3;
d_ext_shaft_disks = [76.2 76.2 50.8]'*1e-3;
m0_disks          = [0 2.54e-4 5.08e-5]./D_disk'*2;
m0_phase          = [0 pi/2 0];

%% Bearing positions
z_bearings       = [81.28 640.08]'*1e-3;
d_extShaft_bear  = [25.4 38.1]'*1e-3;

%% Mesh
[z_coords, idx_sort] = sort([coord_should; z_disks; z_bearings]);
numEl    = length(z_coords)-1;
zero_vec = zeros(numEl,1);

d_extShaft_augm = [d_ExtShaft; d_ext_shaft_disks; d_extShaft_bear];
d_shaftEl = [zero_vec, d_extShaft_augm, zero_vec, d_extShaft_augm];
d_shaftEl = d_shaftEl(idx_sort(2:end)-1,:);

n_intNodes = zero_vec;
n_intNodes([9,12]) = 2; n_intNodes([5,6,7,10,11]) = 3; n_intNodes([1,3,4,8]) = 5;
isGraded          = zero_vec;
shaftElType       = ones(numEl,1);
shaftElProperties = [rho*ones(numEl,1), E*ones(numEl,1), zero_vec];

rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, shaftElType, shaftElProperties);

%% Disk and bearing assembly
rotor.disk(3) = struct('type',[],'node',[],'thick',[],'D_int',[],...
    'D_ext',[],'rho',[],'m0',[],'delta',[]);
[rotor.disk.type] = deal(1); [rotor.disk.rho] = deal(rho);
[rotor.disk.D_int] = num2cell(d_ext_shaft_disks){:};
[rotor.disk.D_ext] = num2cell(D_disk){:};
[rotor.disk.thick] = num2cell(thick_disk){:};
[rotor.disk.m0]    = num2cell(m0_disks){:};
[rotor.disk.delta] = num2cell(m0_phase){:};

z_nodes_mesh = [rotor.nodes.node]; z_coord_mesh = [rotor.nodes.coord];
idx_d = arrayfun(@(z) find(z_coord_mesh==z,1), z_disks);
[rotor.disk.node] = num2cell(idx_d){:};

idx_b1 = find(z_coord_mesh==z_bearings(1),1);
idx_b2 = find(z_coord_mesh==z_bearings(2),1);
k = 5.534e6; c = 1.7513e3;
rotor.bearing(1) = struct('type',3,'node',z_nodes_mesh(idx_b1),'kx',k,'ky',k,'cx',c,'cy',c);
rotor.bearing(2) = struct('type',3,'node',z_nodes_mesh(idx_b2),'kx',k,'ky',k,'cx',c,'cy',c);
rotor.forcing(1).type = 1;

figureRotor(rotor);

%% Campbell diagram
rotorSpeed = (0:1000:6000)*pi/30;
[natPuls, Mode, kappa] = charRoots(rotor, rotorSpeed);
[NatPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);
selectModes = 1:6;
NatPuls_sort0 = NatPuls_sort(selectModes,:);
kappa_sort0   = kappa_sort(:,selectModes,:);
NX = 1; isDamped = true;
plotCampbell(rotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort0);

%% Mode shapes
figure('Name','Modes','NumberTitle','off')
for jj = 1:length(selectModes)
    subplot(2,3,jj); kk = selectModes(jj);
    plotMode(rotor, Mode_sort(:,kk,2), 'number', kk, ...
        'frequency', abs(NatPuls_sort0(1,kk))/2/pi, 'rotorSpeed', rotorSpeed(2)*30/pi);
end

%% Unbalance response
rotorSpeed = (0:10:6000)*pi/30;
response   = FRF(rotor, rotorSpeed);
plotFRF(rotorSpeed*30/pi, response, [24.1 24.2])

%% Critical speed map vs bearing stiffness
stiff = logspace(log10(3e6), log10(2e9), 7);
criticalSpeeds_map = zeros(2, length(stiff));

for jj = 1:length(stiff)
    ks = stiff(jj);
    rotor.bearing(1) = struct('type',3,'node',z_nodes_mesh(idx_b1),'kx',ks,'ky',ks,'cx',c,'cy',c);
    rotor.bearing(2) = struct('type',3,'node',z_nodes_mesh(idx_b2),'kx',ks,'ky',ks,'cx',c,'cy',c);
    [Speeds,~] = critSpeeds(rotor, NX, isDamped, 2, 'num_crit', 4);
    criticalSpeeds_map(:,jj) = Speeds([2,4]);
end

figure('Name','Critical Speed Map','NumberTitle','off','Color','w');
hold on;
loglog(stiff, criticalSpeeds_map(1,:)*30/pi,'b-o','LineWidth',1.5,'MarkerFaceColor','b','MarkerSize',6,'DisplayName','1st Mode');
loglog(stiff, criticalSpeeds_map(2,:)*30/pi,'r-s','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',6,'DisplayName','2nd Mode');
hold off;
set(gca,'XScale','log','YScale','log');
grid on; grid minor;
xlabel('Bearing Stiffness k [N/m]'); ylabel('Critical Speed [rpm]');
title('Critical Speed Map vs Bearing Stiffness'); legend show;

%% Stability analysis with speed-dependent bearings (type 13)
rotor = rmfield(rotor,'bearing');
rotorSpeed_points = [1 1000 3000 5000 7000 10000]*pi/30;
zero_pts = zeros(1,6);
kxx_pts  = [1.58 1.40 1.23 1.05 0.876 0.700]*1e6;
cyy_pts  = [1.75 5.25 8.76 12.3 15.8 19.3]*1e3;

rotor.bearing(1) = struct('type',13,'node',z_nodes_mesh(idx_b1),'Speed_points',rotorSpeed_points,...
    'kxx_points',kxx_pts,'kxy_points',zero_pts,'kyx_points',zero_pts,'kyy_points',kxx_pts,...
    'cxx_points',cyy_pts,'cxy_points',zero_pts,'cyx_points',zero_pts,'cyy_points',cyy_pts);
rotor.bearing(2) = struct('type',13,'node',z_nodes_mesh(idx_b2),'Speed_points',rotorSpeed_points,...
    'kxx_points',kxx_pts,'kxy_points',zero_pts,'kyx_points',zero_pts,'kyy_points',kxx_pts,...
    'cxx_points',cyy_pts,'cxy_points',zero_pts,'cyx_points',zero_pts,'cyy_points',cyy_pts);

RotorSpeed = (10:100:5000)*pi/30;
[natPuls, ~, kappa] = charRoots(rotor, RotorSpeed);
[NatPuls_sort, ~, kappa_sort] = sortModesMAC(natPuls, [], kappa);
NatPuls_sort0 = NatPuls_sort(1:6,:);

for mm = 1:6
    plotRootLocus(RotorSpeed*30/pi, NatPuls_sort0(mm,:));
end
plotCampbell(RotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort(:,1:6,:));

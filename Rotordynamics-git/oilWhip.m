% OILWHIP  Nonlinear analysis of oil whirl and oil whip in a Jeffcott-like
%   rotor on two short fluid-film journal bearings.
%
%   Reproduces a published benchmark case to demonstrate the transition from
%   synchronous response to oil whirl (subsynchronous ~0.5X) and, at higher
%   speeds, oil whip (lock-in to a rotor natural frequency). The rotor is a
%   simple stepped shaft with a single overhung disk and internal material
%   damping (beta = 25e-5 s).
%
%   ROTOR DESCRIPTION
%     Shaft : three sections, total length ≈ 4*Le + 2*L
%             de = 12 mm, Le = 145 mm, L = 10 mm
%     Disk  : D_ext = 94.7 mm, thick = 43 mm, m0 = 75 µg
%             at z = L + 2*Le
%     Bearings: two nonlinear fluid-film (type 7.1) at z = L/2 and
%               z = 1.5*L + 4*Le; D = 20 mm, c = 90 µm, eta = 0.04 Pa·s
%     Forcing : unbalance (type 1) + gravity (type 4)
%     Mesh  : 3 coarse elements refined to 5, 9, 3 sub-elements
%             (cosine grading not used; uniform spacing)
%
%   ANALYSIS PERFORMED
%     1. Rotor schematic (figureRotor)
%     2. Natural frequencies and Campbell diagram at speeds 120–180 Hz
%        with 0.5X and 1X lines (charRoots, sortModesMAC, plotCampbell)
%     3. Root locus (plotRootLocus)
%     4. Speed sweep (180–240 Hz): time simulation at each speed;
%        steady-state response extracted by interpolation onto a uniform
%        time grid (timeSimulation, nr = 4)
%     5. Poincaré map at one speed (poincareMap)
%     6. Bifurcation diagram (plotBifurcation)
%     7. Waterfall (cascade) diagram via DFT of the x-response
%        (DFT, waterfallPlot)
%     8. Results saved to resultsOilWhirl.mat
%
% REFERENCE
%   Based on published experimental/numerical benchmark (see inline comment).
%
% SEE ALSO
%   timeSimulation, nonLinBearingMatrix, plotBifurcation, poincareMap,
%   waterfallPlot, DFT, meshGenerator

clear
close all
addpath(genpath('auxFunc/'))

%% Rotor geometry
de = 0.012;    % shaft diameter [m]
Le = 0.145;    % element length [m]
L  = 0.010;    % bearing span half-length [m]
coord_should = [L; 4*Le+L; 4*Le+2*L];

rho     = 7800;
E       = 206.7e9;
beta    = 25e-5;         % internal damping coefficient [s]
poisson = 0.285;
G       = E / 2 / (1 + poisson);

%% Mesh generation
d_shaftEl     = repmat([0 de 0 de], 3, 1);
z_coords      = [0; coord_should];
n_intNodes    = [1; 3; 1];
isGraded      = zeros(3,1);
shaftElType   = 3*ones(3,1);
shaftElProperties = repmat([rho E poisson], 3, 1);
shaftLoad     = zeros(3,2);
intDamping    = beta * ones(3,1);

rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, ...
    shaftElType, shaftElProperties, shaftLoad, intDamping);

%% Disk
z_disk     = L + 2*Le;
Dext_disk  = 94.7e-3;
thick_disk = 0.043;
m0         = 7.5e-5;  % unbalance mass [kg]

z_coord_mesh = [rotor.nodes.coord];
z_nodes_mesh = [rotor.nodes.node];
idx_disk = find(z_coord_mesh == z_disk, 1);

rotor.disk(1) = struct('type',1,'node',z_nodes_mesh(idx_disk),...
    'D_ext',Dext_disk,'D_int',de,'thick',thick_disk,'rho',rho,'m0',m0);

%% Bearings (nonlinear fluid-film, type 7.1)
rotor.bearing(2) = struct('type',[],'node',[]);
idx_b1 = find(z_coord_mesh == L/2, 1);
idx_b2 = find(z_coord_mesh == 1.5*L + 4*Le, 1);
node_bears_cell = num2cell([z_nodes_mesh(idx_b1), z_nodes_mesh(idx_b2)]);

D_b = L/0.5; c_b = 90e-6; eta_b = 0.04;
[rotor.bearing.type] = deal(7.1);
[rotor.bearing.node] = node_bears_cell{:};
[rotor.bearing.L]    = deal(L);
[rotor.bearing.D]    = deal(D_b);
[rotor.bearing.c]    = deal(c_b);
[rotor.bearing.eta]  = deal(eta_b);

rotor.forcing(1).type = 1;   % mass unbalance
rotor.forcing(2).type = 4;   % gravity

%% Rotor schematic and natural frequencies
figureRotor(rotor);

RotorSpeed = (120:5:180)*2*pi;
[natPuls, Mode, kappa] = charRoots(rotor, RotorSpeed);
[natPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

NX = [0.5 1];
NatPuls_sort0 = natPuls_sort(1:10,:);
kappa_sort0   = kappa_sort(:,1:10,:);
isDamped = true;
plotCampbell(RotorSpeed, NatPuls_sort0, NX, isDamped, kappa_sort0);
plotRootLocus(RotorSpeed*30/pi, NatPuls_sort0(1,:));

%% Speed sweep: time simulations at 180–240 Hz
Rotor_speed = (180:5:240)*2*pi;     % [rad/s]

t_end  = 2*pi/Rotor_speed(1) * 1000;   % 1000 revolutions at lowest speed
t_span = [0, t_end];
t_sampl = linspace(t_end*0.91, t_end, 1000);

niter  = length(Rotor_speed);
numEl_mesh = length(z_nodes_mesh);

resp_x = zeros(1000, niter);
resp_y = zeros(1000, niter);

initPos = zeros(4*numEl_mesh, 1);
initPos(2:4:end) = 1e-6;
initVel = zeros(4*numEl_mesh, 1);

for ii = 1:niter
    [time, q, ~] = timeSimulation(rotor, t_span, Rotor_speed(ii), initPos, initVel, nr=4);

    idx_steady = time > (t_end * 0.9);
    resp_x(:,ii) = interp1(time(idx_steady), q(idx_steady,17), t_sampl, 'linear');
    resp_y(:,ii) = interp1(time(idx_steady), q(idx_steady,18), t_sampl, 'linear');

    fprintf('Step %d/%d completed.\n', ii, niter);
end

%% Post-processing
poincareMap(t_sampl, resp_x(:,11), resp_y(:,11), Rotor_speed(5))
plotBifurcation(Rotor_speed*30/pi, t_sampl', resp_x)

% waterfall diagram
[f, resp_x_freq, ~] = DFT(t_sampl', resp_x, 100000);
idx_freq = f < 400;
epsilon  = 1e-20;
Z_plot   = log10(abs(resp_x_freq(idx_freq,:)) + epsilon)';
waterfallPlot(f(idx_freq), Rotor_speed(:)/(2*pi), Z_plot, [0.5, 1]);

save('resultsOilWhirl.mat')

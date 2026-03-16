clear
close all

addpath(genpath('auxFunc/'))
% recreating example based on 
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
Dext_disk = 94.7e-3; %[m]
thick_disk = 0.043; %[m]
m0 = 7.5e-5; %[Kg]

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

D = L/0.5; %[m]
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
%% natural frequencies
RotorSpeed = (120:5:180)*2*pi;
[natPuls,Mode,kappa] = charRoots(rotor,RotorSpeed);
[natPuls_sort, Mode_sort, kappa_sort] = sortModesMAC(natPuls, Mode, kappa);

%sorted
NX = [0.5 1];
NatPuls_sort0 = natPuls_sort(1:1:10,:);
kappa_sort0 = kappa_sort(:,1:1:10,:);
isDamped = 1;
plotCampbell(RotorSpeed,NatPuls_sort0,NX,isDamped,kappa_sort0);

%% root locus
plotRootLocus(RotorSpeed*30/pi,NatPuls_sort0(1,:));
%% solving for variable clearence
Rotor_speed = (180:5:240)*2*pi ; %[rad/s]

t_end = 2*pi/Rotor_speed(1) * 1000;
t_span = [0,t_end];
t_sampl = linspace(t_end*0.91,t_end,1000);

niter = length(Rotor_speed);
numEl = length(z_nodes_mesh);

resp_x = zeros(1000,niter);
resp_y = zeros(1000,niter);

initPos = zeros(4*(numEl),1);
initPos(2:4:end) = 1e-6;
initVel = zeros(4*(numEl),1);

for ii=1:niter
    [time,q,qdot] = timeSimulation(rotor,t_span,Rotor_speed(ii),initPos,initVel,nr = 4);
  
    idx_steady = time > (t_end * 0.9);
    resp_x(:,ii) = interp1(time(idx_steady),q(idx_steady,17),t_sampl,'linear');
    resp_y(:,ii) = interp1(time(idx_steady),q(idx_steady,18),t_sampl,'linear');

    fprintf('Step %d/%d completato.\n', ii, niter);
end
%% Poincaré map
poincareMap(t_sampl,resp_x(:,11),resp_y(:,11),Rotor_speed(5))

%% Bifurcation map
plotBifurcation(Rotor_speed * 30/pi, t_sampl, resp_x)

%% waterfall plot
[f,resp_x_freq,Nyquist_f] = DFT(t_sampl,resp_x,100000);
idx_freq = f < 400;
raw_Z = abs(resp_x_freq(idx_freq,:));
epsilon = 1e-20; 

Z_plot = log10(raw_Z + epsilon)'; % Matrice [Speed x Freq]
X_freq = f(idx_freq);
Y_speed_Hz = Rotor_speed(:) / (2*pi);
NX = [0.5,1];

waterfallPlot(X_freq,Y_speed_Hz,Z_plot,NX);
save('resultsOilWhirl.mat')
function [] = plotMode(Rotor,Mode,options)
%PLOTMODE plots the Mode Shape for a set of Modes from charRoots
%
%INPUT: Rotor       Structure
%       Mode        Array of eigevectors
%
%   plotMode(rotor,Mode,'number',[],'frequency',[],'rotorSpeed',[]);
%
%   number = number associated to the mode
%   max_iter = maximum interations number
%   eigenvalue = frequency associated to the mode [Hz]
arguments (Input)
    Rotor (1,1) struct
    Mode (:,:) double {mustBeVector}
    options.number (1,1) double {mustBeInteger, mustBePositive} = NaN
    options.frequency (1,1) double {mustBePositive} = NaN
    options.rotorSpeed (1,1) double = NaN
end
nodes = Rotor.nodes;
n_nodes = numel(nodes);

if isscalar([nodes.node])
    error('Mode cannot be plotted')
end
j = 1i;
ndof = 4 * n_nodes;

modex = Mode(1:4:ndof);
modey = Mode(2:4:ndof);
z = [Rotor.nodes.coord];

% define circles
np = 200;
circle = exp(j*(0:np)*2*pi/np).';

x_matrix = real(modex * circle.');
y_matrix = real(modey * circle.');

xp_matrix = [zeros(1, n_nodes); x_matrix.'; zeros(1, n_nodes)];
yp_matrix = [zeros(1, n_nodes); y_matrix.'; zeros(1, n_nodes)];
zp_matrix = ones(np + 3, 1) * z;

xp = xp_matrix(:);
yp = yp_matrix(:);
zp = zp_matrix(:);

% define the line using shape functions
nn = 20; 
zeta = (0:nn).'/nn;
onn = ones(nn+1,1);
N1 = onn - 3*zeta.^2 + 2*zeta.^3;
N2 = zeta - 2*zeta.^2 + zeta.^3;
N3 = 3*zeta.^2 - 2*zeta.^3;
N4 = -zeta.^2 + zeta.^3;

n_points_line = (nn + 1) * (n_nodes - 1);
xn = zeros(n_points_line, 1);
yn = zeros(n_points_line, 1);
zn = zeros(n_points_line, 1);

for inode = 1:n_nodes-1
   idx_start = (inode - 1) * (nn + 1) + 1;
   idx_end = inode * (nn + 1);
   indices = idx_start:idx_end;
   
   Le = z(inode+1)-z(inode);
   
   xx = [N1, Le*N2, N3, Le*N4] * real(Mode([4*inode-3 4*inode 4*inode+1 4*inode+4]));
   yy = [N1, -Le*N2, N3, -Le*N4] * real(Mode([4*inode-2 4*inode-1 4*inode+2 4*inode+3],1));   
   zz = z(inode)*onn+Le*zeta;
   
   xn(indices) = xx;
   yn(indices) = yy;
   zn(indices) = zz;
end

max_xy = max(abs([xp; yp]));
if max_xy == 0
    max_xy = 1; 
end
xn = xn/max_xy;
yn = yn/max_xy;
xp = xp/max_xy;
yp = yp/max_xy;

max_z = max(abs(zn));
if max_z == 0
    max_z = 1;
end
zn = zn/max_z;
zp = zp/max_z;

hl = plot3(zp,xp,yp,'-k',zn,xn,yn,'-k');
set( hl(1), 'LineWidth', 0.5)
set( hl(2), 'LineWidth', 1.5)

sf = 2;

xlim([0 1]);
ylim([-sf sf]);
zlim([-sf sf]);

pbaspect([1.5 1 1]);
set(gca, 'Color', 'none');

freq = options.frequency;
number = options.number;
rotorSpeed = options.rotorSpeed;

line1 = []; line2 = []; line3=[];
if ~isnan(number),line1 = sprintf('Mode: %d ', number);end
if ~isnan(rotorSpeed),line2 = sprintf('Frequency: %.2f [Hz]',freq);end
if ~isnan(freq), line3 = sprintf('Rotor Speed: %.2f [rpm]', rotorSpeed);end
sgtitle('Non-dimensional Modes Representation')
title({line1, line2, line3});
end
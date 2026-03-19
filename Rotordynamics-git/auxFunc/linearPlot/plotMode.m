function [] = plotMode(Rotor, Mode, options)
% PLOTMODE  Plots a non-dimensional 3D representation of a rotor mode shape.
%
%   For each node, the precession orbit is drawn as a circle in the (x, y)
%   plane at the corresponding axial position z, parameterised by the
%   complex amplitudes of the mode shape. The orbit is defined by the
%   complex exponential sweep:
%       x(t) = Re(modex * exp(i*alpha)),  alpha in [0, 2*pi)
%       y(t) = Re(modey * exp(i*alpha))
%
%   The deflected shaft centreline is also drawn using Hermite cubic shape
%   functions (real part of the mode shape). All coordinates are
%   normalised to unit maximum displacement for shape comparison.
%
% SYNTAX
%   plotMode(Rotor, Mode)
%   plotMode(Rotor, Mode, "number", n, "frequency", f, "rotorSpeed", Omega)
%
% INPUT ARGUMENTS
%   Rotor - (1x1 struct) Rotor data structure with at least:
%             .nodes  - node array with .node and .coord [m]
%   Mode  - (ndof x 1 complex double) Mode shape eigenvector as returned
%           by charRoots (physical partition, ndof entries).
%           DOF ordering: [u1, v1, theta_u1, theta_v1, u2, ...].
%           The function uses x-DOFs (indices 1:4:ndof) and y-DOFs
%           (indices 2:4:ndof) for the orbit, and shape functions for the
%           shaft centreline.
%
% NAME-VALUE OPTIONS
%   "number"      - (positive integer) Mode number for the title annotation
%   "frequency"   - (positive double) Natural frequency [Hz] for annotation
%   "rotorSpeed"  - (double) Rotor speed [rpm] for the title annotation
%
% OUTPUT
%   A 3D figure is added to the current figure (or a new one if none
%   exists) showing:
%     - Thin black lines : nodal precession orbits (circles at each node)
%     - Thick black line : deflected shaft centreline (Hermite interpolation)
%   The z-axis represents the normalised axial coordinate, the x- and
%   y-axes the normalised lateral displacements. An optional title reports
%   mode number, frequency, and rotor speed if provided.
%
% NOTES
%   - Requires at least two nodes; raises an error for single-node rotors.
%   - The mode shape is normalised by the maximum nodal orbit amplitude so
%     that the plot is always unit-scale, regardless of eigenvector scaling.
%   - This function does not create a new figure; call figure() beforehand
%     if needed (e.g., for a subplot layout).
%
% EXAMPLE
%   speeds = linspace(0, 1500, 100);
%   [evals, evecs] = charRoots(Rotor, speeds);
%   figure;
%   plotMode(Rotor, evecs(:, 1), "number", 1, ...
%            "frequency", abs(imag(evals(1,50)))/(2*pi), ...
%            "rotorSpeed", speeds(50)*30/pi);
%
% SEE ALSO
%   charRoots, plotOrbit, plotDisplacement, computeDisplacement

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

% nodal orbit circles (parameterised by alpha in [0, 2*pi])
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

% shaft centreline via Hermite cubic shape functions
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

% normalise all coordinates to unit amplitude
max_xy = max(abs([xp; yp]));
if max_xy == 0
    max_xy = 1; 
end
xn = xn/max_xy; yn = yn/max_xy;
xp = xp/max_xy; yp = yp/max_xy;

max_z = max(abs(zn));
if max_z == 0
    max_z = 1;
end
zn = zn/max_z;
zp = zp/max_z;

hl = plot3(zp,xp,yp,'-k',zn,xn,yn,'-k');
set( hl(1), 'LineWidth', 0.5)   % orbit circles
set( hl(2), 'LineWidth', 1.5)   % shaft centreline

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
if ~isnan(number),     line1 = sprintf('Mode: %d ', number); end
if ~isnan(rotorSpeed), line2 = sprintf('Frequency: %.2f [Hz]', freq); end
if ~isnan(freq),       line3 = sprintf('Rotor Speed: %.2f [rpm]', rotorSpeed); end
sgtitle('Non-dimensional Modes Representation')
title({line1, line2, line3});
end

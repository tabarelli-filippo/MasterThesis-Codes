function [] = plotStressAnalysis(Rotor, response, nu)
% PLOTSTRESSANALYSIS  Visualises the results of a Von Mises static stress
%   analysis along the deflected rotor shaft as a 3D colour-mapped plot.
%
%   The deflected shaft centreline is reconstructed from the nodal
%   displacement vector using Hermite cubic shape functions (20
%   intermediate points per element). The static safety factor nu is
%   linearly interpolated along the shaft axis and mapped onto the
%   deflected line using a colour scale. Node positions and their safety
%   factor values are annotated on the plot.
%
%   The colour scale uses an inverted jet colormap clipped to [1.5, 10]:
%     - Red  (nu ≈ 1.5) : low safety margin — near yielding
%     - Blue (nu ≈ 10)  : high safety margin
%
% SYNTAX
%   plotStressAnalysis(Rotor, response, nu)
%
% INPUT ARGUMENTS
%   Rotor    - (1x1 struct) Rotor data structure with at least:
%                .nodes  - node array with .node and .coord [m]
%                .shaft  - shaft element array (geometry fields used for
%                          displacement reconstruction via shape functions)
%   response - (ndof x 1 double or complex double) Nodal displacement
%              vector. Only the real part is used. Typically the output of
%              FRF at a specific speed or a static deflection vector.
%              DOF ordering: [u1, v1, theta_u1, theta_v1, ...]
%   nu       - ((n_nodes) x 1 double) Static safety factor at each node,
%              as returned by stressAnalysis. nu = sigma_y / sigma_VM.
%
% OUTPUT
%   A 3D figure titled 'Stress Analysis - Von Mises Criterion' is created
%   showing:
%     - Colour-mapped thick line : deflected shaft centreline, coloured
%       by the local safety factor nu
%     - Dash-dot line  : undeflected (reference) shaft axis
%     - Black circles  : nodal positions annotated with nu values
%   The colorbar label indicates the static safety coefficient nu_st.
%   Axes: z [m] (axial), x [mm] and y [mm] (lateral, converted for
%   visibility; displacement values are scaled ×1000 from [m] to [mm]).
%
% EXAMPLE
%   resp = FRF(Rotor, omega_crit);
%   [nu_s, ~, ~] = stressAnalysis(Rotor, real(resp), 250e6);
%   plotStressAnalysis(Rotor, resp, nu_s);
%
% SEE ALSO
%   stressAnalysis, FRF, plotDisplacement, computeDisplacement

arguments (Input)
    Rotor (1,1) struct
    response (:,:) {mustBeVector} 
    nu (:,:) double {mustBeVector}
end

n_nodes = numel(Rotor.nodes);
nodes = 1:1:n_nodes;
z = [Rotor.nodes.coord];

% Hermite cubic shape functions
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
   
   Le = z(inode+1) - z(inode);
   
   q_u = real(response([4*inode-3, 4*inode, 4*inode+1, 4*inode+4]));
   q_v = real(response([4*inode-2, 4*inode-1, 4*inode+2, 4*inode+3]));
   
   xx_m = [N1, Le*N2, N3, Le*N4] * q_u;
   yy_m = [N1, -Le*N2, N3, -Le*N4] * q_v; 
   
   xn(indices) = xx_m * 1000;   % convert [m] to [mm]
   yn(indices) = yy_m * 1000;
   
   zn(indices) = z(inode)*onn + Le*zeta;
end

% interpolate safety factor along shaft axis
nu_interp = interp1(z, nu, zn, 'linear', 'extrap');

% nodal positions along interpolated curve
idx_nodes = (nodes - ones(size(nodes))) * (nn + 1); idx_nodes(1) = 1;
x_nodes = xn(idx_nodes);
y_nodes = yn(idx_nodes);

%% plot
figure('Name', 'Rotor Stress Analysis','NumberTitle', 'off');
hold on
% deflected shaft centreline, colour-mapped by safety factor
surface([zn,zn], [xn,xn], [yn,yn], [nu_interp,nu_interp], ...
    'FaceColor', 'none','EdgeColor', 'interp','LineWidth', 4);
% undeformed centreline
plot3([zn(1), zn(end)], [0, 0], [0, 0] , 'k-.', 'LineWidth', 0.5, ...
    'DisplayName', 'Undeformed Centerline');

% node markers and safety factor annotations
plot3(z', x_nodes, y_nodes,'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
nu_text = compose('%.2f', nu);
text(z, x_nodes, y_nodes, nu_text,'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom','FontSize', 10);

view(3);
grid on;
grid minor;
colormap(flipud(jet));
c = colorbar;
c.Label.String = 'Static safety coefficient \nu_{st}';
clim([1.5,10]);

xlabel('z [m]'); ylabel('x [mm]'); zlabel('y [mm]');
pbaspect([1 0.25 0.25]);
set(gca, 'Color', 'none');
title('Stress Analysis - Von Mises Criterion');
hold off
end

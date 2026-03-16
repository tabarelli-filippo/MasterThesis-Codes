function [] = plotStressAnalysis(Rotor,response,nu)
%PLOTSTRESSANALYSIS Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    Rotor (1,1) struct
    response (:,:) {mustBeVector} 
    nu (:,:) double {mustBeVector}
end

n_nodes = numel(Rotor.nodes);
nodes = 1:1:n_nodes;
z = [Rotor.nodes.coord];

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
   
   Le = z(inode+1) - z(inode); % Element length
   
   q_u = real(response([4*inode-3, 4*inode, 4*inode+1, 4*inode+4]));
   q_v = real(response([4*inode-2, 4*inode-1, 4*inode+2, 4*inode+3]));
   
   xx_m = [N1, Le*N2, N3, Le*N4] * q_u;
   yy_m = [N1, -Le*N2, N3, -Le*N4] * q_v; 
   
   xn(indices) = xx_m * 1000 ;
   yn(indices) = yy_m * 1000 ;
   
   zn(indices) = z(inode)*onn + Le*zeta;
end

nu_interp = interp1(z,nu,zn,'linear', 'extrap');
% nodes coordinates
idx_nodes = (nodes - ones(size(nodes))) * (nn + 1); idx_nodes(1) = 1;
x_nodes = xn(idx_nodes);
y_nodes = yn(idx_nodes);

%% plot
figure('Name', 'Rotor Stress Analysis','NumberTitle', 'off');
hold on
% rotor line
surface([zn,zn], [xn,xn], [yn,yn], [nu_interp,nu_interp],'FaceColor', 'none','EdgeColor', 'interp','LineWidth', 4);
% undeformed centerline
plot3([zn(1), zn(end)], [0, 0], [0, 0] , 'k-.', 'LineWidth', 0.5, 'DisplayName', 'Undeformed Centerline');

% nodes and values
plot3(z', x_nodes, y_nodes,'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
nu_text = compose('%.2f', nu);
text(z, x_nodes, y_nodes, nu_text,'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','FontSize', 10);

view(3);
grid on;
grid minor;
colormap(flipud(jet));
c = colorbar;
c.Label.String = 'Static security coefficient \nu_{st}';
clim([1.5,10]);

xlabel('z [m]'); ylabel('x [mm]'); zlabel('y [mm]');
pbaspect([1 0.25 0.25]);
set(gca, 'Color', 'none');
title('Stress Analysis - Von Mises Criterion');
hold off
end
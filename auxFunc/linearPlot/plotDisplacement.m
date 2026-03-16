function [] = plotDisplacement(Rotor, Disp,options)
%PLOTDISPLACEMENT plots displacement in both xz and yx planes. Also maximum
%   displacement is plotted
%
%INPUT: Rotor       structure of the rotor
%       Disp        complex response to evaluate
%       rotorSpeed  rotor speed for plotting text    
arguments (Input)
    Rotor (1,1) struct
    Disp (:,:) double {mustBeVector}
    options.rotorSpeed (1,1) double {mustBePositive} = nan
end

nodes = Rotor.nodes;
n_nodes = numel(nodes);

if isscalar([nodes.node])
    error('Displacement cannot be plotted')
end
ndof = 4 * n_nodes;
z = [Rotor.nodes.coord];

% shape functions
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
   
   xx = [N1, Le*N2, N3, Le*N4] * real(Disp([4*inode-3 4*inode 4*inode+1 4*inode+4]));
   yy = [N1, -Le*N2, N3, -Le*N4] * real(Disp([4*inode-2 4*inode-1 4*inode+2 4*inode+3],1));   
   zz = z(inode)*onn+Le*zeta;
   
   xn(indices) = xx;
   yn(indices) = yy;
   zn(indices) = zz;
end

%% nodes
idx_x_nodes = 1:4:ndof;
idx_y_nodes = 2:4:ndof;

x_nodes = real(Disp(idx_x_nodes));
y_nodes = real(Disp(idx_y_nodes));
z_nodes = z(:); % Ensure column vector

maxDisp_smooth = sqrt(xn.^2 + yn.^2); 
maxDisp_nodes  = sqrt(x_nodes.^2 + y_nodes.^2); 

%% plot
[global_max_val, max_idx] = max(maxDisp_smooth);
z_at_max = zn(max_idx); 
limit_val = max(abs([xn(:); yn(:); maxDisp_smooth(:)])) * 1.1;

figure('Name', 'Displacement Profile', 'NumberTitle', 'off', 'Color', 'w');
t = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

if ~isnan(options.rotorSpeed)
    titleStr = sprintf('Lateral Displacement - Rotor speed: %.2f [rpm]', options.rotorSpeed);
else
    titleStr = 'Rotordynamic Lateral Displacement';
end
title(t, titleStr, 'FontSize', 14, 'FontWeight', 'bold');

markerStyle = {'ko', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', 'MarkerSize', 5, 'LineStyle', 'none'};
zeroLineStyle = {'-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0};
zeros_at_nodes = zeros(size(z_nodes));

% x displacement
nexttile
hold on;
yline(0, zeroLineStyle{:}); 
plot(z_nodes, zeros_at_nodes, markerStyle{:});
plot(zn, xn, 'k-', 'LineWidth', 1.2); 
plot(z_nodes, x_nodes, markerStyle{:});        

grid on; grid minor;
xlabel('Axial Position [m]');
ylabel('x-Disp [m]');
ylim([-limit_val, limit_val]);

% y displacement
nexttile
hold on;
yline(0, zeroLineStyle{:});
plot(z_nodes, zeros_at_nodes, markerStyle{:});
plot(zn, yn, 'k-', 'LineWidth', 1.2); 
plot(z_nodes, y_nodes, markerStyle{:});        

grid on; grid minor;
xlabel('Axial Position [m]');
ylabel('y-Disp [m]');
ylim([-limit_val, limit_val]);


% Magnitude
nexttile
hold on;
yline(0, zeroLineStyle{:});
plot(z_nodes, zeros_at_nodes, markerStyle{:});
plot(zn, maxDisp_smooth, 'k-', 'LineWidth', 1.5); 
plot(z_nodes, maxDisp_nodes, markerStyle{:});              

% Highlight Max Value
plot(z_at_max, global_max_val, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
text_str = sprintf(' Max: %.2e', global_max_val);
text(z_at_max, global_max_val, text_str, 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left', 'FontSize', 9, 'FontWeight','bold', ...
    'Color', 'r');

grid on; grid minor;
xlabel('Axial Position [m]');
ylabel('|Disp| [m]');
ylim([0, limit_val]);


linkaxes(findall(gcf, 'type', 'axes'), 'x');
xlim([min(zn) max(zn)]);
end
function [] = figureRotor(Rotor)
% FIGUREROTOR  Generates a 2D schematic cross-section of the rotor system,
%   showing the shaft profile, disks, and bearing positions.
%
%   The figure provides a rapid visual check of the rotor geometry as
%   defined in the Rotor structure. Shaft sections are drawn as filled
%   outlines (external and internal profiles), disks are shown as
%   trapezoidal cross-sections, and bearings are represented as filled
%   triangular symbols above and below the shaft centre-line.
%   Shaft steps (diameter discontinuities at shared nodes) are drawn
%   explicitly as vertical lines.
%
% SYNTAX
%   figureRotor(Rotor)
%
% INPUT ARGUMENTS
%   Rotor - (1x1 struct) Rotor data structure. Required fields:
%             .nodes   - array of node structs with .node and .coord [m]
%             .shaft   - array of shaft element structs with:
%                          .d_ext1, .d_ext2  external diameters [m]
%                          .d_int1, .d_int2  internal diameters [m]
%             .disk    - array of disk structs; type-dependent fields:
%                          Type 1: .D_ext, .D_int, .thick [m]
%                          Type 2: .Id, .mass (used to estimate radius)
%             .bearing - array of bearing structs with .node and .type
%
% OUTPUT
%   A figure titled 'Rotor Schematic' is created with:
%     - Black outline : shaft external and internal profiles [mm]
%     - Red lines     : disk cross-sections
%     - Blue triangles: bearing symbols (above and below shaft axis)
%     - Dash-dot line : shaft rotation axis
%     - Node numbers  annotated at each node
%     - Bearing type  annotated below each bearing symbol
%
% NOTES
%   - All dimensions are converted to millimetres for display.
%   - For disk type 2, the external radius is estimated from the diametral
%     moment of inertia: R_ext ≈ sqrt(4*Id/m).
%   - Hollow shaft internal profiles are shown as dashed lines only when
%     the internal diameter is non-zero.
%
% EXAMPLE
%   figureRotor(Rotor);
%
% SEE ALSO
%   plotMode, plotDisplacement, plotStressAnalysis, meshGenerator

arguments (Input)
    Rotor (1,1) struct
end
%% line widths
disk_line = 1;
draw_line = 0.75;
axis_line = 1.5;

%% structures
Nodes = Rotor.nodes; Shaft = Rotor.shaft; Disk = Rotor.disk; Bearing= Rotor.bearing;

n_nodes = numel(Nodes); n_shaft = numel(Shaft); n_disks = numel(Disk); n_bearing = numel(Bearing);

figure('Name','Rotor Schematic','NumberTitle', 'off')
hold on

%% plot nodes
nodes_coord = [Nodes.coord]*1e3;
z_coord = sort([nodes_coord, nodes_coord(2:end-1)]);

plot(z_coord,zeros(size(z_coord)),'ro')
% node labels
for ii = 1:n_nodes
    text(nodes_coord(ii), 0, num2str(Nodes(ii).node),'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','FontSize', 8);
end
%% plot shaft
d_ext1 = [Shaft.d_ext1]*1e3;
d_ext2 = [Shaft.d_ext2]*1e3;
if isempty([Shaft.d_ext2])
    d_ext2 = d_ext1;
end
d_ext = zeros(size(z_coord));

d_ext(1:2:end) = d_ext1/2;
d_ext(2:2:end) = d_ext2/2;

d_int1 = [Shaft.d_int1]*1e3;
d_int2 = [Shaft.d_int2]*1e3;

if isempty([Shaft.d_int2])
    d_int2 = d_int1;
end
d_int = nan(size(z_coord));

d_int(1:2:end) = d_int1/2;
d_int(2:2:end) = d_int2/2;

z_poly = [0,z_coord, fliplr(z_coord),0];
r_poly_e = [0,d_ext, -fliplr(d_ext),0];
r_poly_i = [d_int, -fliplr(d_int)];

% external shaft profile
plot(z_poly, r_poly_e, 'k', 'LineWidth', draw_line);
% internal (hollow) profile, shown only when present
if ~isequal(d_int1,zeros(size(d_int1)))
    plot(z_poly(2:end-1), r_poly_i, 'k--', 'LineWidth', 0.5);
end

% rotation axis
plot([-r_poly_e(2), z_coord(end)+r_poly_e(2)], [0, 0], 'k-.', 'LineWidth',axis_line);
% shaft diameter steps (shoulders)
for i = 1:1:(length(z_coord) - 1)
    if (d_ext(i) ~= d_ext(i+1)) && (z_coord(i) == z_coord(i+1))
        r_min = min(d_ext(i), d_ext(i+1));    
        plot([z_coord(i), z_coord(i)], [-r_min, r_min], 'k','LineWidth', draw_line); 
    end
end
%% bearings
r_shaft_aux = [d_ext1,0;0,d_ext2]/2;
bearing_nodes = [Bearing.node];
bearing_coord = nodes_coord(bearing_nodes);
h_bearings = zeros(1,length(bearing_nodes));

for ii = 1:length(bearing_nodes)
    x_bear = bearing_coord(ii);
    y_bear = max(r_shaft_aux(:,bearing_nodes(ii)));
    aa = y_bear/4;
    x1 = x_bear - aa; x2 = x_bear + aa; y1 = y_bear + aa;
    h_bearings(ii)= fill([x_bear x1 x2 x_bear],[y_bear y1 y1 y_bear],'b','EdgeColor','b','LineWidth',draw_line);
    fill([x_bear x1 x2 x_bear],-1*[y_bear y1 y1 y_bear],'b','EdgeColor','b','LineWidth',draw_line)
    bear_type = sprintf('Bearing Type: %.1f',Bearing(ii).type);
    text(x_bear, -y1, bear_type,'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','FontSize', 8);
end

%% plot disks
disk_nodes = [Disk.node];
z_disks = nodes_coord(disk_nodes);
type = [Disk.type];
h_disks = zeros(1,length(disk_nodes));
R_ext_mem = zeros(1,length(disk_nodes));
thick_mem = zeros(1,length(disk_nodes));

for ii = 1:length(disk_nodes)
    switch type(ii)
        case 1
            R_ext = Disk(ii).D_ext*0.5e3;
            R_int = Disk(ii).D_int*0.5e3;
            half_thick = Disk(ii).thick*0.5e3;
        case 2
            R_ext = sqrt(4*Disk(ii).Id/Disk(ii).mass)*1e3;
            R_int = max(r_shaft_aux(:,disk_nodes(ii)));
            half_thick = 0.2*R_ext;
    end
    
    x1 = z_disks(ii)-half_thick;
    x2 = z_disks(ii)+half_thick;
    y1 = R_int + half_thick*tan(10*pi/180);
    if R_ext < y1
        R_ext = R_ext + y1;
    end

    R_ext_mem(ii) = R_ext;
    thick_mem(ii) = half_thick;

    x_disk = [z_disks(ii) x1 x1 x2 x2 z_disks(ii) ];
    y_disk = [R_int, y1 , R_ext, R_ext, y1 , R_int];
    h_disks(ii) = plot(x_disk,y_disk,'r','LineWidth',disk_line);
    plot(x_disk,-y_disk,'r','LineWidth',disk_line)
end

%% finishing
sgtitle('ROTOR SCHEMATIC');

xlabel('Axial Distance [mm]');
ylabel('Shaft Radius [mm]');
y_max = max(R_ext_mem*1.25);
x_min = max(r_poly_e(2),thick_mem(1));
x_max = max(r_poly_e(2),thick_mem(end));
xlim([-x_min*1.1,z_coord(end)+x_max*1.1]);
ylim([-y_max,y_max]);

grid on; grid minor;

handle = [h_bearings(1), h_disks(1)];
leg_text = {'Bearings','Disks'};
legend(handle,leg_text);
set(gca, 'Color', 'none');

legend('show')
hold off
end

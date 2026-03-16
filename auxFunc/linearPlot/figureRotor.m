function [] = figureRotor(Rotor)
%FIGUREROTOR plots a schematic of the rotor system
%
%INPUT: Rotor   Structure of the rotor
arguments (Input)
    Rotor (1,1) struct
end
%% lines
disk_line = 1;
draw_line = 0.75;
axis_line = 1.5;

%% structures
Nodes = Rotor.nodes; Shaft = Rotor.shaft; Disk = Rotor.disk; Bearing= Rotor.bearing;

n_nodes = numel(Nodes); n_shaft = numel(Shaft); n_disks = numel(Disk); n_bearing = numel(Bearing);

figure('Name','Rotor Schematic','NumberTitle', 'off')
hold on

%% plot nodes
% extract nodes
nodes_coord = [Nodes.coord]*1e3;
z_coord = sort([nodes_coord, nodes_coord(2:end-1)]);

plot(z_coord,zeros(size(z_coord)),'ro')
% Node id
for ii = 1:n_nodes
    text(nodes_coord(ii), 0, num2str(Nodes(ii).node),'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','FontSize', 8);
end
%% plot shaft 
%extract external diameters
d_ext1 = [Shaft.d_ext1]*1e3;
d_ext2 = [Shaft.d_ext2]*1e3;
if isempty([Shaft.d_ext2])
    d_ext2 = d_ext1;
end
d_ext = zeros(size(z_coord));

d_ext(1:2:end) = d_ext1/2;
d_ext(2:2:end) = d_ext2/2;

%extract internal diameters
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

% shaft
plot(z_poly, r_poly_e, 'k', 'LineWidth', draw_line);
if ~isequal(d_int1,zeros(size(d_int1)))
    plot(z_poly(2:end-1), r_poly_i, 'k--', 'LineWidth', 0.5);
end

% dash dotted line
plot([-r_poly_e(2), z_coord(end)+r_poly_e(2)], [0, 0], 'k-.', 'LineWidth',axis_line);
%plot shaft shoulders
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

    R_ext_mem(ii) = R_ext; %for graph scaling
    thick_mem(ii) = half_thick; %for graph scaling

    x_disk = [z_disks(ii) x1 x1 x2 x2 z_disks(ii) ];
    y_disk = [R_int, y1 , R_ext, R_ext, y1 , R_int];
    h_disks(ii) = plot(x_disk,y_disk,'r','LineWidth',disk_line);
    plot(x_disk,-y_disk,'r','LineWidth',disk_line)
end



%% finishing
sgtitle('ROTOR SCHEMATIC');

xlabel('Axial Distance [mm]');
ylabel('Shaft Radius [mm]');
y_max = max(R_ext_mem*1.25); % Calculate maximum height for y-axis limits
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
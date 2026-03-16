function [Rotor] = meshGenerator(z_coords,d_shaftEl,n_intNodes,isGraded,shaftElType,...
    shaftElProperties,shaftLoad,intDamping)
%MESHGENERATOR creates the structure with mesh definition. Each element can
%be refined with mode points as defined.

arguments (Input)
    z_coords (:,1) double
    d_shaftEl (:,4) double
    n_intNodes (:,1) double
    isGraded (:,1) double
    shaftElType (:,1) double
    shaftElProperties (:,3) double
    shaftLoad(:,2) double = nan(1,1)
    intDamping (:,1) double = nan(1,1)
end


numEl_original = length(z_coords) - 1;
numNodi_original = length(z_coords);

if ~isnan(shaftLoad)
    axForce = shaftLoad(:,1);
    torque = shaftLoad(:,2);
else
    axForce = zeros(numEl_original,1);
    torque = zeros(numEl_original,1);
end

if isnan(intDamping)
    intDamping = zeros(numEl_original,1);
end

% elements per node
num_sub_elements = n_intNodes + 1;
idx_rep = repelem((1:numEl_original)', num_sub_elements);

shaftElType_mesh = shaftElType(idx_rep, :);
shaftElProperties_mesh = shaftElProperties(idx_rep, :);
axForce_mesh = axForce(idx_rep, :);
torque_mesh = torque(idx_rep, :);
intDamping_mesh = intDamping(idx_rep, :);

%% initialization
total_nodes = numNodi_original + sum(n_intNodes);
z_coords_mesh = zeros(total_nodes, 1);
z_coords_mesh(1) = z_coords(1);

total_new_elements = sum(num_sub_elements);
P = size(d_shaftEl, 2);
P_half = P / 2;
d_shaftEl_mesh = zeros(total_new_elements, P);

current_node_idx = 1;
current_elem_idx = 0; 
%% adding new coordinates
for ii = 1:numEl_original
    n1_z = z_coords(ii);
    n2_z = z_coords(ii+1);
    
    nn = n_intNodes(ii);
    n_sub_el = num_sub_elements(ii);
    n_points = nn + 2;
    
    if isGraded(ii) % graded distribution
        points = ((1 - cos(linspace(0, pi, n_points))) / 2)';
    else            % linear distribution
        points = linspace(0, 1, n_points)';
    end
    
    coord_new_nodes = n1_z + points(2:end) * (n2_z - n1_z);
    
    start_idx = current_node_idx + 1;
    end_idx = current_node_idx + n_sub_el;
    z_coords_mesh(start_idx:end_idx) = coord_new_nodes;
    current_node_idx = end_idx;
    % interpolates diameters
    props_start_orig = d_shaftEl(ii, 1:P_half);
    props_end_orig = d_shaftEl(ii, (P_half + 1):P);
    % interpolate properties over adimesional points
    interp_props_at_nodes = interp1([0; 1], [props_start_orig; props_end_orig], points, 'linear');
    % re-define shaftEl
    sub_elements_props = [interp_props_at_nodes(1:end-1, :), interp_props_at_nodes(2:end, :)];
    
    start_idx = current_elem_idx + 1;
    end_idx = current_elem_idx + n_sub_el;
    % refined mesh
    d_shaftEl_mesh(start_idx:end_idx, :) = sub_elements_props;
    current_elem_idx = end_idx;
end

numEl = length(z_coords_mesh);
%% Assembla l'output
rho = shaftElProperties_mesh(:,1); 
E = shaftElProperties_mesh(:,2); 
poisson = shaftElProperties_mesh(:,3); 
G = E./2./(1 + poisson);

%% Diameter definitions
dint1 = d_shaftEl_mesh(:,1); dext1 = d_shaftEl_mesh(:,2); 
dint2 = d_shaftEl_mesh(:,3); dext2 = d_shaftEl_mesh(:,4);

%%  Rotor definition
Rotor.nodes(numEl) = struct('node',[],'coord',[]);
Rotor.shaft(numEl-1) = struct('type', [], 'node1', [], 'node2', [], 'd_int1', ...
    [], 'd_ext1', [],'d_int2',[],'d_ext2',[],'rho', [], 'E', [], 'G', []);

%% Rotor.nodes
nnodes = 1:numEl;
nodes_cell = num2cell(nnodes);
[Rotor.nodes.node] = nodes_cell{:};
z_nodes_cell = num2cell(z_coords_mesh);
[Rotor.nodes.coord] = z_nodes_cell{:};

%% Rotor.shaft
shaftElType_cell = num2cell(shaftElType_mesh);
[Rotor.shaft.type] = shaftElType_cell{:};

rho_cell = num2cell(rho);
[Rotor.shaft.rho] = rho_cell{:};

E_cell = num2cell(E);
[Rotor.shaft.E] = E_cell{:};

G_cell = num2cell(G);
[Rotor.shaft.G] = G_cell{:};

%nodes
[Rotor.shaft.node1] = nodes_cell{1:end-1};
[Rotor.shaft.node2] = nodes_cell{2:end};

%external diameter
Dint1_cell = num2cell(dint1);
[Rotor.shaft.d_int1] = Dint1_cell{:};

Dext1_cell = num2cell(dext1);
[Rotor.shaft.d_ext1] = Dext1_cell{:};

Dint2_cell = num2cell(dint2);
[Rotor.shaft.d_int2] = Dint2_cell{:};

Dext2_cell = num2cell(dext2);
[Rotor.shaft.d_ext2] = Dext2_cell{:};

%% load definitions
axForce_cell = num2cell(axForce_mesh);
[Rotor.shaft.AxialForce] = axForce_cell{:};

torque_cell = num2cell(torque_mesh);
[Rotor.shaft.Torque] = torque_cell{:};

%% internal damping definition
intDamping_cell = num2cell(intDamping_mesh);
[Rotor.shaft.beta] = intDamping_cell{:};
end

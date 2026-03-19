function [Rotor] = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, shaftElType, ...
    shaftElProperties, shaftLoad, intDamping)
% MESHGENERATOR  Builds the rotor finite element mesh (Rotor structure)
%   from a coarse geometric description, with optional element refinement.
%
%   Each coarse element (defined by two consecutive z_coords entries) is
%   subdivided into (n_intNodes + 1) sub-elements. Node coordinates can be
%   distributed uniformly or with a cosine grading (denser near element
%   ends) to better capture stress concentrations and mode shape curvature.
%   Shaft cross-section diameters are linearly interpolated at the new
%   intermediate nodes to preserve the taper profile.
%
% SYNTAX
%   Rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, ...
%               shaftElType, shaftElProperties)
%   Rotor = meshGenerator(z_coords, d_shaftEl, n_intNodes, isGraded, ...
%               shaftElType, shaftElProperties, shaftLoad, intDamping)
%
% INPUT ARGUMENTS
%   z_coords          - (N+1 x 1 double) Axial coordinates of the coarse
%                       mesh nodes [m]. N = number of coarse elements.
%                       Must be monotonically increasing.
%   d_shaftEl         - (N x 4 double) Shaft diameters at each coarse
%                       element in the order:
%                         [d_int1, d_ext1, d_int2, d_ext2]
%                       where subscripts 1 and 2 denote node 1 and node 2
%                       of each element [m]. Use d_int = 0 for solid shaft.
%   n_intNodes        - (N x 1 double) Number of intermediate nodes to
%                       insert in each coarse element. n_intNodes(i) = 0
%                       means no refinement (1 sub-element = original).
%   isGraded          - (N x 1 logical/double) Refinement type per element:
%                         0 → uniform spacing
%                         1 → cosine grading (denser at ends)
%   shaftElType       - (N x 1 double) Element type for each coarse element
%                       (1 = Euler-Bernoulli, 2/3 = Timoshenko; see
%                       shaftElement for the full list). Replicated to
%                       all sub-elements of each coarse element.
%   shaftElProperties - (N x 3 double) Material properties per coarse
%                       element: [rho [kg/m³], E [Pa], Poisson_ratio [-]].
%                       G is computed as E / (2*(1+nu)).
%   shaftLoad         - (N x 2 double, optional, default NaN) Mechanical
%                       loads per coarse element: [AxialForce [N], Torque [N·m]].
%                       Pass nan(1,1) or omit to apply no loads.
%   intDamping        - (N x 1 double, optional, default NaN) Internal
%                       damping coefficient beta [s] per coarse element
%                       (proportional damping: C_int = beta * K).
%                       Pass nan(1,1) or omit for undamped elements.
%
% OUTPUT ARGUMENTS
%   Rotor - (struct) Rotor data structure with fields:
%             .nodes  - array of node structs (.node, .coord)
%             .shaft  - array of refined shaft element structs
%                       (.type, .node1, .node2, .d_ext1, .d_int1,
%                        .d_ext2, .d_int2, .rho, .E, .G,
%                        .AxialForce, .Torque, .beta)
%           Note: .disk and .bearing fields are NOT populated by this
%           function and must be added manually after calling meshGenerator.
%
% NOTES
%   - Cosine grading uses the formula: zeta = (1 - cos(linspace(0,pi,n+2)))/2,
%     which produces denser node spacing near element endpoints.
%   - Diameter values are linearly interpolated at sub-element nodes.
%   - G is derived from E and Poisson's ratio: G = E / (2*(1+nu)).
%   - The output Rotor.shaft arrays are of length (total sub-elements).
%
% EXAMPLE
%   % Two-element shaft: uniform cylindrical + tapered section
%   z   = [0; 0.3; 0.7];              % coarse node coordinates [m]
%   d   = [0, 0.06, 0, 0.06;         % element 1: solid, uniform d=60 mm
%          0, 0.06, 0, 0.04];         % element 2: solid, tapered 60→40 mm
%   n_int = [3; 5];                   % 3 and 5 intermediate nodes
%   grad  = [0; 1];                   % uniform; cosine
%   type  = [2; 21];                  % Timoshenko; tapered Euler-Bernoulli
%   prop  = [7800, 2.1e11, 0.3;       % steel
%            7800, 2.1e11, 0.3];
%   Rotor = meshGenerator(z, d, n_int, grad, type, prop);
%   % Then add disks and bearings:
%   Rotor.disk(1) = struct('type', 2, 'node', 4, 'mass', 5, ...);
%   Rotor.bearing(1) = struct('type', 1, 'node', 1);
%
% SEE ALSO
%   rotorMatrix, shaftElement, shaftTaperEl, figureRotor, massRotor

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

numEl_original   = length(z_coords) - 1;
numNodi_original = length(z_coords);

if ~isnan(shaftLoad)
    axForce = shaftLoad(:,1);
    torque  = shaftLoad(:,2);
else
    axForce = zeros(numEl_original,1);
    torque  = zeros(numEl_original,1);
end

if isnan(intDamping)
    intDamping = zeros(numEl_original,1);
end

% replicate coarse-element properties to all sub-elements
num_sub_elements = n_intNodes + 1;
idx_rep = repelem((1:numEl_original)', num_sub_elements);

shaftElType_mesh       = shaftElType(idx_rep, :);
shaftElProperties_mesh = shaftElProperties(idx_rep, :);
axForce_mesh           = axForce(idx_rep, :);
torque_mesh            = torque(idx_rep, :);
intDamping_mesh        = intDamping(idx_rep, :);

%% Initialise refined mesh arrays
total_nodes        = numNodi_original + sum(n_intNodes);
z_coords_mesh      = zeros(total_nodes, 1);
z_coords_mesh(1)   = z_coords(1);

total_new_elements = sum(num_sub_elements);
P      = size(d_shaftEl, 2);
P_half = P / 2;
d_shaftEl_mesh = zeros(total_new_elements, P);

current_node_idx = 1;
current_elem_idx = 0; 

%% Insert intermediate nodes and interpolate diameters
for ii = 1:numEl_original
    n1_z = z_coords(ii);
    n2_z = z_coords(ii+1);
    
    nn      = n_intNodes(ii);
    n_sub_el = num_sub_elements(ii);
    n_points = nn + 2;
    
    if isGraded(ii)   % cosine grading: denser at element ends
        points = ((1 - cos(linspace(0, pi, n_points))) / 2)';
    else              % uniform spacing
        points = linspace(0, 1, n_points)';
    end
    
    coord_new_nodes = n1_z + points(2:end) * (n2_z - n1_z);
    
    start_idx = current_node_idx + 1;
    end_idx   = current_node_idx + n_sub_el;
    z_coords_mesh(start_idx:end_idx) = coord_new_nodes;
    current_node_idx = end_idx;

    % linearly interpolate diameters at sub-node positions
    props_start_orig = d_shaftEl(ii, 1:P_half);
    props_end_orig   = d_shaftEl(ii, (P_half + 1):P);
    interp_props_at_nodes = interp1([0; 1], [props_start_orig; props_end_orig], points, 'linear');

    % re-pack as [d_start, d_end] per sub-element
    sub_elements_props = [interp_props_at_nodes(1:end-1, :), interp_props_at_nodes(2:end, :)];
    
    start_idx = current_elem_idx + 1;
    end_idx   = current_elem_idx + n_sub_el;
    d_shaftEl_mesh(start_idx:end_idx, :) = sub_elements_props;
    current_elem_idx = end_idx;
end

numEl = length(z_coords_mesh);

%% Compute material properties from input
rho     = shaftElProperties_mesh(:,1); 
E       = shaftElProperties_mesh(:,2); 
poisson = shaftElProperties_mesh(:,3); 
G       = E ./ 2 ./ (1 + poisson);

%% Diameter definitions (node 1 and node 2 of each sub-element)
dint1 = d_shaftEl_mesh(:,1); dext1 = d_shaftEl_mesh(:,2); 
dint2 = d_shaftEl_mesh(:,3); dext2 = d_shaftEl_mesh(:,4);

%%  Assemble Rotor structure
Rotor.nodes(numEl) = struct('node',[],'coord',[]);
Rotor.shaft(numEl-1) = struct('type', [], 'node1', [], 'node2', [], ...
    'd_int1', [], 'd_ext1', [], 'd_int2', [], 'd_ext2', [], ...
    'rho', [], 'E', [], 'G', []);

%% Rotor.nodes
nnodes     = 1:numEl;
nodes_cell = num2cell(nnodes);
[Rotor.nodes.node]  = nodes_cell{:};
z_nodes_cell = num2cell(z_coords_mesh);
[Rotor.nodes.coord] = z_nodes_cell{:};

%% Rotor.shaft — element type and material
shaftElType_cell = num2cell(shaftElType_mesh);
[Rotor.shaft.type] = shaftElType_cell{:};

rho_cell = num2cell(rho);
[Rotor.shaft.rho] = rho_cell{:};

E_cell = num2cell(E);
[Rotor.shaft.E] = E_cell{:};

G_cell = num2cell(G);
[Rotor.shaft.G] = G_cell{:};

% connectivity
[Rotor.shaft.node1] = nodes_cell{1:end-1};
[Rotor.shaft.node2] = nodes_cell{2:end};

% diameters
Dint1_cell = num2cell(dint1); [Rotor.shaft.d_int1] = Dint1_cell{:};
Dext1_cell = num2cell(dext1); [Rotor.shaft.d_ext1] = Dext1_cell{:};
Dint2_cell = num2cell(dint2); [Rotor.shaft.d_int2] = Dint2_cell{:};
Dext2_cell = num2cell(dext2); [Rotor.shaft.d_ext2] = Dext2_cell{:};

%% Rotor.shaft — loads and internal damping
axForce_cell = num2cell(axForce_mesh);
[Rotor.shaft.AxialForce] = axForce_cell{:};

torque_cell = num2cell(torque_mesh);
[Rotor.shaft.Torque] = torque_cell{:};

intDamping_cell = num2cell(intDamping_mesh);
[Rotor.shaft.beta] = intDamping_cell{:};
end

function [xn,yn,zn] = computeDisplacement(Rotor,Disp,options)
% COMPUTEDISPLACEMENT  Reconstructs the continuous lateral displacement
%   field of the rotor shaft from the nodal displacement vector using
%   Hermite cubic shape functions.
%
%   Intermediate points are evaluated within each shaft element by
%   interpolating the four nodal DOFs (displacement and rotation at each
%   end node) with the standard Euler-Bernoulli beam shape functions.
%   This allows smooth visualisation of the deflected shaft shape.
%
% SYNTAX
%   [xn, yn, zn] = computeDisplacement(Rotor, Disp)
%   [xn, yn, zn] = computeDisplacement(Rotor, Disp, nn=N)
%
% INPUT ARGUMENTS
%   Rotor - (1x1 struct) Rotor data structure containing at least:
%             .nodes  - array of node structs with fields .node and .coord
%             .shaft  - shaft element definitions (used implicitly via nodes)
%   Disp  - (ndof x 1 complex double) Nodal displacement vector. Only
%           the real part is used for plotting. Typically the output of
%           FRF or a modal eigenvector.
%           DOF ordering: [u1, v1, theta_u1, theta_v1, u2, ...] where
%           u is the x-displacement and v is the y-displacement.
%   nn    - (name-value, positive integer, default 20) Number of
%           intermediate evaluation points per element. Higher values
%           yield smoother shaft-line plots.
%
% OUTPUT ARGUMENTS
%   xn - ((nn+1)*(n_nodes-1) x 1 double) x-displacement at interpolated
%        points along the shaft axis [m]
%   yn - ((nn+1)*(n_nodes-1) x 1 double) y-displacement at interpolated
%        points along the shaft axis [m]
%   zn - ((nn+1)*(n_nodes-1) x 1 double) Axial coordinate of each
%        interpolated point along the shaft [m]
%
% SHAPE FUNCTIONS
%   The four Hermite cubic shape functions are parameterised by the
%   normalised coordinate zeta in [0, 1] over element length Le:
%     N1(zeta) = 1 - 3*zeta^2 + 2*zeta^3
%     N2(zeta) = zeta - 2*zeta^2 + zeta^3      (multiplied by Le)
%     N3(zeta) = 3*zeta^2 - 2*zeta^3
%     N4(zeta) = -zeta^2 + zeta^3              (multiplied by Le)
%
% NOTES
%   - The function requires at least two nodes; it raises an error for
%     single-node rotors.
%   - Only the real part of Disp is used (steady-state amplitude plot).
%   - For orbit plots (whirl visualisation), call this function separately
%     for the real and imaginary parts of the complex response.
%
% EXAMPLE
%   resp = FRF(Rotor, omega_crit);        % complex response at critical speed
%   [xn, yn, zn] = computeDisplacement(Rotor, resp, nn=50);
%   figure;
%   plot3(zn, xn, yn); grid on;
%   xlabel('z [m]'); ylabel('x [m]'); zlabel('y [m]');
%   title('Deflected shaft shape at critical speed');
%
% SEE ALSO
%   FRF, charRoots, rotorMatrix

arguments (Input)
    Rotor (1,1) struct
    Disp (:,:) double {mustBeVector}
    options.nn (1,1) double {mustBeInteger, mustBePositive} = 20
end

nn = options.nn;
nodes = Rotor.nodes;
n_nodes = numel(nodes);

if isscalar([nodes.node])
    error('Displacement cannot be plotted')
end

z = [Rotor.nodes.coord];

% shape functions
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

end

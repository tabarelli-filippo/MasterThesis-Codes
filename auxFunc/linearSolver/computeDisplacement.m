function [xn,yn,zn] = computeDisplacement(Rotor,Disp,options)
%COMPUTEDISPLACEMENT: computes displacement in both planes at each node.
%Intermediate points can be evaluated.
%
%INPUT: Rotor   structure of the rotor
%       Disp    complex response to evaluate
%       nn      number of intermediate points for shape functions
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
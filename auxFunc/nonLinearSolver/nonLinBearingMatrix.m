function [Force] = nonLinBearingMatrix(Rotor,rotorSpeed,q,qdot)
%BEARINGMATRIX  calculates the mass, stiffness and damping matrices for the
% bearings. Also gives the list of constrained dofs
%
% INPUT:Rotor       Structure
%       rotorSpeed  Double ONLY
%       q           position vector
%       qdot        velocity vector
%
% OUTPUT:   Mb  Bearing Mass matrix
%           Cb  Bearing Damping matrix
%           Kb  Bearing Stiffness matrix
%           zero_dof  Constrained Dofs

nodes = Rotor.nodes;
bearings = Rotor.bearing;

% initialize matrices
n_nodes  = numel(nodes);
nBearing = numel(bearings);
ndof = 4 * n_nodes;


Force = zeros(ndof,1);

for ii = 1:nBearing
    type = bearings(ii).type;
    nnode = bearings(ii).node;
    dof = (4*nnode-3):4*nnode;
    u = real(q(4*nnode-3));
    v = real(q(4*nnode-2));
    udot = real(qdot(4*nnode-3));
    vdot = real(qdot(4*nnode-2));
    bearForce = zeros(4,1);
    switch type
        case 7.1 % nonLinear Oil Bearing
            if isfield(Rotor.bearing(ii), 'D'), D = Rotor.bearing(ii).D; end % bearing diameter [m]
            if isfield(Rotor.bearing(ii), 'L'), L = Rotor.bearing(ii).L; end % bearing length [m]
            if isfield(Rotor.bearing(ii), 'c'), c = Rotor.bearing(ii).c; end % bearing radial clearance [m]
            if isfield(Rotor.bearing(ii), 'eta'), eta = Rotor.bearing(ii).eta; end % oil viscosity [Ns/m^2]
            
            if rotorSpeed == 0
                error('>>>> Error - bearing type 7 - fluid bearing model undefined at zero speed')
            end
            
            aa = v*rotorSpeed + 2*udot;
            bb = u*rotorSpeed - 2*vdot;
            cuv= c^2-u^2-v^2;

            alpha = atan(aa/bb) - pi/2*sign(aa/bb) - pi/2*sign(aa);
            G = 2*c/sqrt(cuv)*(pi/2 + atan((v*cos(alpha)-u*sin(alpha))/sqrt(cuv)));
            S = c*(u*cos(alpha) + v*sin(alpha))/(c^2 - (u*cos(alpha) + v*sin(alpha))^2);
            V = (2*c^2 + c*(v*cos(alpha) - u*sin(alpha))*G)/cuv;

            dd1 = eta/8*rotorSpeed*L^(3)*D/c^2;
            dd2 = -c*sqrt(bb^2+aa^2)/rotorSpeed/cuv;

            fxx = dd1*dd2*(3*V*u/c - G*sin(alpha) - 2*S*cos(alpha));
            fyy = dd1*dd2*(3*V*v/c + G*cos(alpha) - 2*S*sin(alpha));
            
            bearForce = [fxx;fyy;0;0];
            
        case 8.1 %Non Linear Seals
    end
    Force(dof) = Force(dof) + bearForce;
end
end
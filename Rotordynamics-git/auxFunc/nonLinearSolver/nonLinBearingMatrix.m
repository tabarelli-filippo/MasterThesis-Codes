function [Force] = nonLinBearingMatrix(Rotor, rotorSpeed, q, qdot)
% NONLINBEARINGMATRIX  Computes the nonlinear bearing reaction force vector
%   for use in time-domain rotor simulations.
%
%   Unlike bearingMatrix (which returns linearised K and C matrices valid
%   for frequency-domain analysis), this function evaluates the actual
%   nonlinear bearing force at each time step given the instantaneous shaft
%   position and velocity. It is called internally by the ODE right-hand
%   side functions in timeSimulation and runUp.
%
%   Currently implemented nonlinear bearing types:
%     Type 7.1 - Nonlinear short journal bearing (Reynolds equation,
%                Ocvirk short-bearing theory). The hydrodynamic film
%                force is a nonlinear function of the instantaneous
%                eccentricity (u, v) and velocity (udot, vdot).
%     Type 8.1 - Nonlinear seal (reserved; not yet implemented)
%
% SYNTAX
%   Force = nonLinBearingMatrix(Rotor, rotorSpeed, q, qdot)
%
% INPUT ARGUMENTS
%   Rotor      - (struct) Rotor data structure. Bearing fields required
%                for type 7.1:
%                  .bearing(i).D   bearing diameter [m]
%                  .bearing(i).L   bearing length [m]
%                  .bearing(i).c   radial clearance [m]
%                  .bearing(i).eta oil dynamic viscosity [N·s/m²]
%   rotorSpeed - (scalar double) Instantaneous shaft angular velocity
%                [rad/s]. Must be > 0 (fluid film undefined at rest).
%   q          - (ndof x 1 double) Instantaneous nodal displacement
%                vector [m]. Only translational DOFs u and v at the
%                bearing node are used. DOF ordering:
%                [u1, v1, theta_u1, theta_v1, u2, ...]
%   qdot       - (ndof x 1 double) Instantaneous nodal velocity vector
%                [m/s]. Same DOF ordering as q.
%
% OUTPUT ARGUMENTS
%   Force - (ndof x 1 double) Global nonlinear bearing force vector [N].
%           Non-zero only at the translational DOFs (u, v) of bearing
%           nodes. Rotational DOFs receive zero force contribution.
%
% THEORY (Type 7.1)
%   The short-bearing nonlinear force components are derived from the
%   Reynolds equation under the Ocvirk (short-bearing) approximation.
%   The film pressure is integrated analytically over the half-film
%   (pi-film cavitation model), yielding closed-form expressions for
%   fx and fy as functions of (u, v, udot, vdot, Omega). The attitude
%   angle alpha is computed from the squeeze and wedge velocity components.
%
% NOTES
%   - For linear (small-displacement) analysis, use bearingMatrix instead.
%   - rotorSpeed = 0 causes a division by zero and triggers an error.
%   - ndof = 4 * n_nodes
%
% EXAMPLE
%   % Typically called internally by runUp or timeSimulation:
%   F_bearing = nonLinBearingMatrix(Rotor, Omega, q_current, qdot_current);
%
% SEE ALSO
%   bearingMatrix, timeSimulation, runUp

nodes = Rotor.nodes;
bearings = Rotor.bearing;

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
        case 7.1 % nonlinear short journal bearing (Ocvirk pi-film model)
            if isfield(Rotor.bearing(ii), 'D'), D = Rotor.bearing(ii).D; end
            if isfield(Rotor.bearing(ii), 'L'), L = Rotor.bearing(ii).L; end
            if isfield(Rotor.bearing(ii), 'c'), c = Rotor.bearing(ii).c; end
            if isfield(Rotor.bearing(ii), 'eta'), eta = Rotor.bearing(ii).eta; end
            
            if rotorSpeed == 0
                error('>>>> Error - bearing type 7.1 - fluid bearing model undefined at zero speed')
            end
            
            % squeeze and wedge velocity components
            aa = v*rotorSpeed + 2*udot;
            bb = u*rotorSpeed - 2*vdot;
            cuv = c^2 - u^2 - v^2;

            % attitude angle and auxiliary integrals (Ocvirk solution)
            alpha = atan(aa/bb) - pi/2*sign(aa/bb) - pi/2*sign(aa);
            G = 2*c/sqrt(cuv)*(pi/2 + atan((v*cos(alpha)-u*sin(alpha))/sqrt(cuv)));
            S = c*(u*cos(alpha) + v*sin(alpha))/(c^2 - (u*cos(alpha) + v*sin(alpha))^2);
            V = (2*c^2 + c*(v*cos(alpha) - u*sin(alpha))*G)/cuv;

            dd1 = eta/8*rotorSpeed*L^3*D/c^2;
            dd2 = -c*sqrt(bb^2+aa^2)/rotorSpeed/cuv;

            fxx = dd1*dd2*(3*V*u/c - G*sin(alpha) - 2*S*cos(alpha));
            fyy = dd1*dd2*(3*V*v/c + G*cos(alpha) - 2*S*sin(alpha));
            
            bearForce = [fxx; fyy; 0; 0];
            
        case 8.1 % nonlinear seal — reserved for future implementation
    end
    Force(dof) = Force(dof) + bearForce;
end
end

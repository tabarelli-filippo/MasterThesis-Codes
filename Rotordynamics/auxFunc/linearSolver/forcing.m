function [unForce,bendForce,gravityAcc] = forcing(Rotor)
% FORCING  Assembles the complex forcing vectors for steady-state response
%   analysis of a rotor-bearing system.
%
%   The total forcing in the frequency domain is:
%       F(Omega) = bendForce + unForce * Omega^2 + gravityAcc
%
%   Four excitation types are supported and can be combined:
%
%   Type 1 - Mass unbalance: centrifugal force due to a residual mass
%     eccentricity on one or more disks. The force amplitude scales with
%     Omega^2 and is applied to the translational DOFs of the disk node.
%     Disk fields: .m0 (unbalance mass [kg]) and .D_ext (diameter [m]),
%     OR .epsilon (eccentricity [m], requires the disk mass from diskElement).
%     Optional: .delta (phase angle [rad], default 0).
%
%   Type 2 - Slope (angular) unbalance: moment due to a tilt of the disk
%     principal axis with respect to the shaft axis. Applied to the
%     rotational DOFs. Force amplitude scales with Omega^2.
%     Disk fields: .beta (tilt angle [rad]), .gamma (phase [rad]).
%
%   Type 3 - Shaft bow (static bend): equivalent static force vector
%     induced by a permanent shaft bow. The bow is defined by prescribing
%     displacements at selected master nodes; the slave DOF forces are
%     computed via static condensation of the shaft stiffness matrix.
%     Rotor fields: .bend(i).node, .bend(i).x, .bend(i).y.
%
%   Type 4 - Gravity: constant gravitational body force applied in the
%     negative y-direction at all translational DOFs (g = 9.80665 m/s²).
%
% SYNTAX
%   [unForce, bendForce, gravityAcc] = forcing(Rotor)
%
% INPUT ARGUMENTS
%   Rotor - (1x1 struct) Rotor data structure. Must contain:
%             .nodes    - node definitions
%             .forcing  - array of structs with field .type (1, 2, 3, or 4)
%             .disk     - disk definitions (required for types 1 and 2)
%             .bend     - bow definitions (required for type 3)
%
% OUTPUT ARGUMENTS
%   unForce    - (ndof x 1 complex double) Speed-independent part of the
%                unbalance forcing. Multiply by Omega^2 to get the actual
%                force at speed Omega [N or N·m]
%   bendForce  - (ndof x 1 complex double) Equivalent nodal force vector
%                for the shaft bow excitation [N or N·m]
%   gravityAcc - (ndof x 1 double) Gravity acceleration vector [m/s²]
%                (non-zero only at translational y-DOFs)
%
% NOTES
%   - All force vectors are complex to represent amplitude and phase.
%   - DOF ordering: [u1, v1, theta_u1, theta_v1, u2, ...] where
%     u is the x-direction and v is the y-direction displacement.
%   - ndof = 4 * n_nodes
%   - Multiple forcing types can coexist in Rotor.forcing simultaneously.
%   - The gravity vector is constant (not speed-dependent) and is handled
%     separately in the response solver (FRF).
%
% EXAMPLE
%   [uF, bF, gA] = forcing(Rotor);
%   % Evaluate total force at 1000 rad/s
%   Omega = 1000;
%   F_total = bF + uF * Omega^2;
%
% SEE ALSO
%   FRF, diskElement, rotorMatrix

arguments (Input)
    Rotor (1,1) struct
end
j = 1i;
n_nodes = numel(Rotor.nodes);
ndof = 4 * n_nodes;

[~,~,~,K0,~] = rotorMatrix(Rotor);

forceType = [Rotor.forcing.type];
nforce = length(forceType);
unForce = zeros(ndof,1);
bendForce = zeros(ndof,1);
gravityAcc = zeros(ndof,1);

for nn = 1:nforce
    switch forceType(nn)
        case 1
            %% mass unbalancing
            Disk = Rotor.disk;
            n_disks = numel(Disk);

            for ii = 1:n_disks
                %extracting values
                node = Disk(ii).node;

                delta = 0;
                unbal_mag = 0;

                if isfield(Disk(ii), 'm0')
                    m0 = Disk(ii).m0;
                    D_ext = Disk(ii).D_ext;
                    unbal_mag = 0.5*m0*D_ext;
                end

                if isfield(Disk(ii), 'epsilon')
                    epsilon = Disk(ii).epsilon;
                    [M0e,~] = diskElement(Disk(ii));
                    m = M0e(1,1);
                    unbal_mag = m*epsilon;
                end

                if isfield(Disk(ii), 'delta')
                    delta = Disk(ii).delta;
                end

                force_dof = [4*node-3; 4*node-2];

                unForce(force_dof) = unForce(force_dof) + unbal_mag*exp(j*delta)*[1; -j];
            end
        case 2
            %% slope unbalancing
            Disk = Rotor.disk;
            n_disks = numel(Disk);

            for ii = 1:n_disks
                %extracting values
                node = Disk(ii).node;

                gamma = 0;
                beta = 0;

                if isfield(Disk(ii), 'beta')
                    beta = Disk(ii).beta;
                end

                if isfield(Disk(ii), 'gamma')
                    gamma = Disk(ii).gamma;
                end

                [M0e,C1e] = diskElement(Disk(ii));
                Id = M0e(4,4); Ip = C1e(3,4);
                force_dof = [4*node-1; 4*node];
                unbal_mag = beta*(Id - Ip);
                unForce(force_dof) = unForce(force_dof) + unbal_mag*exp(j*gamma)*[j; 1];
            end
        case 3
            %% Shaft bending (bow)
            if ~isfield(Rotor,'bend')
                error('Bending forcing - Bending not defined')
            end
            n_bend = numel(Rotor.bend); % number of bending points defined

            bendNodes = [Rotor.bend.node];
            bend_x = [Rotor.bend.x];
            bend_y = [Rotor.bend.y];
            bendDef = [bendNodes' bend_x' bend_y'];

            masterdof = zeros(1, 2 * n_bend);
            xmaster = complex(zeros(2 * n_bend, 1));
            for ibend = 1:n_bend
                inode = bendDef(ibend,1);

                idx = (ibend-1)*2 + 1 : ibend*2;

                masterdof(idx) = [4*inode-3, 4*inode-2];
                xnode = bendDef(ibend,2)*[1; -1i] + bendDef(ibend,3)*[1i; 1];
                xmaster(idx) = xnode;
            end
            slavedof = 1:ndof; slavedof(masterdof) = [];
            Kss = K0(slavedof,slavedof);
            Ksm = K0(slavedof,masterdof);
            xslave = - Kss\(Ksm*xmaster);
            xbend = zeros(ndof,1);
            xbend(slavedof) = xslave;
            xbend(masterdof) = xmaster;
            bendForce = K0*xbend;
        case 4
            %% Gravity
            g = -9.80665; %[m/s^2]
            gravityAcc(2:4:end) = g;
    end
end
end

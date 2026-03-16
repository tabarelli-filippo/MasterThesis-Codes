function [unForce,bendForce,gravityAcc] = forcing(Rotor)
%FORCING computes forcing vector from definition, vector is complex
% force = bendForce + unForce*rotorSpeed^2
%
%INPUT:   Rotor       Structure
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
            %% Shaft bending
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
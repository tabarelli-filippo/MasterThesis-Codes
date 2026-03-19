function [Mb,Cb,Kb,zero_dof, eccentricity] = bearingMatrix(Rotor,rotorSpeed)
% BEARINGMATRIX  Assembles the global bearing mass, damping, and stiffness
%   matrices for all bearings defined in the rotor structure. Also returns
%   the list of kinematically constrained (zeroed) degrees of freedom and,
%   for fluid-film bearings, the eccentricity ratio at the given speed.
%
%   The function supports thirteen bearing/support types:
%     Type  1 - Pinned support (constrains translational DOFs)
%     Type  2 - Clamped support (constrains all 4 DOFs at the node)
%     Type  3 - Linear bearing: diagonal K and C, translational DOFs only
%     Type  4 - Linear bearing: diagonal K and C, all 4 DOFs
%     Type  5 - Linear bearing: full 2x2 K and C matrices, no rotations
%     Type  6 - Linear bearing: full 4x4 K and C matrices
%     Type  7 - Short fluid-film journal bearing (Ocvirk's short-bearing
%               theory); K and C computed from the eccentricity ratio
%     Type  8 - Annular pressure seal (Childs' bulk-flow model)
%     Type  9 - Ball/rolling-element bearing, simplified 2-DOF model
%     Type 10 - Squeeze-film damper with optional centering spring
%     Type 11 - Thomas-Alford cross-coupling force (turbomachinery)
%     Type 12 - Labyrinth seal (Scharrer-Childs model)
%     Type 13 - Speed-dependent bearing: K and C interpolated from tables
%
% SYNTAX
%   [Mb, Cb, Kb, zero_dof, eccentricity] = bearingMatrix(Rotor, rotorSpeed)
%
% INPUT ARGUMENTS
%   Rotor       - (struct) Rotor data structure. Must contain the fields:
%                   .nodes    - array of node structs with coordinates
%                   .bearing  - array of bearing structs; each must have:
%                                 .type  (integer, see list above)
%                                 .node  (node index)
%                                 plus type-specific fields (see below)
%   rotorSpeed  - (double, scalar) Shaft rotational speed [rad/s].
%                 Required for speed-dependent bearing types (7, 8, 10, 13).
%
% OUTPUT ARGUMENTS
%   Mb          - (ndof x ndof double) Global bearing mass matrix
%   Cb          - (ndof x ndof double) Global bearing damping matrix
%   Kb          - (ndof x ndof double) Global bearing stiffness matrix
%   zero_dof    - (1 x nzero double) Indices of constrained DOFs
%                 (set to zero displacement by bearing types 1 and 2)
%   eccentricity- (nBearing x 1 double) Eccentricity ratio for each
%                 bearing; non-zero only for fluid-film bearings (type 7)
%
% BEARING-SPECIFIC FIELDS
%   Type 3  : .kx, .ky, .cx, .cy  (scalars, [N/m] and [N·s/m])
%   Type 4  : .k, .c              (4x1 vectors)
%   Type 5  : .K, .C              (2x2 matrices)
%   Type 6  : .K, .C              (4x4 matrices)
%   Type 7  : .F   static load [N]
%             .D   bearing diameter [m]
%             .L   bearing length [m]
%             .c   radial clearance [m]
%             .eta oil dynamic viscosity [N·s/m²]
%             .Linear (optional, default 0) set to 1 to suppress K and C
%   Type 8  : .P   pressure drop across seal [Pa]
%             .R   seal radius [m]
%             .L   seal length [m]
%             .c   radial clearance [m]
%             .V   mean axial flow velocity [m/s]
%             .fric friction coefficient [-]
%   Type 9  : .n_ball  number of balls (8, 12 or 16)
%             .d_ball  ball diameter [m]
%             .fs    static load [N]
%             .alpha contact angle [rad]
%   Type 10 : .eta  oil viscosity [N·s/m²]
%             .R    damper radius [m]
%             .L    damper length [m]
%             .c    clearance [m]
%             .ks   (optional) centering spring stiffness [N/m]
%   Type 11 : .beta empirical cross-coupling coefficient [-]
%             .T    stage torque [N·m]
%             .D    blade mean diameter [m]
%             .L    blade radial length [m]
%   Type 12 : .sealType, .P_reserv, .P_sump, .Nt, .Cr, .L, .T, .Rs,
%             .B, .U_inlet, .omega, .nu  (see labySeals for details)
%   Type 13 : .Speed_points  speed vector [rad/s]
%             .kxx_points, .kxy_points, .kyx_points, .kyy_points
%             .cxx_points, .cxy_points, .cyx_points, .cyy_points
%
% NOTE
%   The global DOF ordering is [u, v, theta_u, theta_v] at each node,
%   yielding ndof = 4 * n_nodes total DOFs.
%
% EXAMPLE
%   [Mb, Cb, Kb, zero_dof, ecc] = bearingMatrix(Rotor, 1000);
%
% SEE ALSO
%   rotorMatrix, charRoots, FRF, critSpeeds

nodes = Rotor.nodes;
bearings = Rotor.bearing;

% initialize matrices
n_nodes  = numel(nodes);
nBearing = numel(bearings);
ndof = 4 * n_nodes;


Mb = zeros(ndof,ndof);
Cb = zeros(ndof,ndof);
Kb = zeros(ndof,ndof);
zero_dof = [];

eccentricity = zeros(nBearing,1);

for ii = 1:nBearing
    type = bearings(ii).type;

    Kb1 = zeros(4,4);
    Cb1 = zeros(4,4);
    Mb1 = zeros(4,4);
    switch type
        case 1 % pinned bearing
            n1 = bearings(ii).node;
            zero_dof = [zero_dof 4*n1-3 4*n1-2];

        case 2 % long, stiff bearing - clamped boundary condition
            n1 = bearings(ii).node;
            zero_dof = [zero_dof 4*n1-3:4*n1];

        case 3       % constant stiffness and damping, diagonal, no rotations
            
            kx=0; ky=0; cx=0; cy=0;

            if isfield(Rotor.bearing(ii), 'kx'), kx = Rotor.bearing(ii).kx; end
            if isfield(Rotor.bearing(ii), 'ky'), ky = Rotor.bearing(ii).ky; end
            if isfield(Rotor.bearing(ii), 'cx'), cx = Rotor.bearing(ii).cx; end
            if isfield(Rotor.bearing(ii), 'cy'), cy = Rotor.bearing(ii).cy; end

            Kb1 = diag( [kx ky 0 0] );
            Cb1 = diag( [cx cy 0 0] );

        case 4       % constant stiffness and damping, diagonal
            if isfield(Rotor.bearing(ii), 'k')
                k = Rotor.bearing(ii).k;
            else
                k = zeros(4,1);
            end

            if isfield(Rotor.bearing(ii), 'c')
                c = Rotor.bearing(ii).c;
            else
                c = zeros(4,1);
            end

            Kb1 = diag(k);
            Cb1 = diag(c);

        case 5       % constant stiffness and damping, no rotations
            if isfield(Rotor.bearing(ii), 'K')
                K = Rotor.bearing(ii).K;
            else
                K = zeros(2,2);
            end

            if isfield(Rotor.bearing(ii), 'C')
                C = Rotor.bearing(ii).C;
            else 
                C = zeros(2,2);
            end

            Z = zeros(2,2);
            Kb1 = [ K, Z ; Z , Z ];
            Cb1 = [ C, Z ; Z , Z ];

        case 6      % constant stiffness and damping, full 4x4 matrices required
            if isfield(Rotor.bearing(ii), 'K')
                K = Rotor.bearing(ii).K;
            else
                K = zeros(4,4);
            end
            if isfield(Rotor.bearing(ii), 'C')
                C = Rotor.bearing(ii).C;
            else
                C = zeros(4,4);
            end

            Kb1 = K;
            Cb1 = C;

        case 7      % fluid film bearings - Linear behaviour
            if isfield(Rotor.bearing(ii), 'F'), F = Rotor.bearing(ii).F; end % static load [N]
            if isfield(Rotor.bearing(ii), 'D'), D = Rotor.bearing(ii).D; end % bearing diameter [m]
            if isfield(Rotor.bearing(ii), 'L'), L = Rotor.bearing(ii).L; end % bearing length [m]
            if isfield(Rotor.bearing(ii), 'c'), c = Rotor.bearing(ii).c; end % bearing radial clearance [m]
            if isfield(Rotor.bearing(ii), 'eta'), eta = Rotor.bearing(ii).eta; end % oil viscosity [Ns/m^2]
            
            Linear = 0;
            if isfield(Rotor.bearing(ii), 'Linear'), Linear = Rotor.bearing(ii).Linear; end % 1 if not linear

            if rotorSpeed == 0
                error('>>>> Error - bearing type 7 - fluid bearing model undefined at zero speed')
            end

            % Find roots of quartic in n^2  (n is eccentricity ratio)
            H = ( 8*c^2*F/(D*rotorSpeed*eta*L^3) )^2;
            n2all = sort( roots([1 -4 (6-(16-pi^2)/H) -(4+pi^2/H) 1]) );
            nroot = 0; n2 = []; % test roots - eccentricity should be between 0 and 1

            % check if solution exists
            for ir=1:4
                nn = n2all(ir);
                if nn>0 && nn<1 && isreal(nn)
                    nroot = nroot + 1;
                    n2 = [n2; nn];
                end
            end

            if nroot == 0
                n2 = 0.5;
                disp('Error in calculating fluid bearing coefficents - no solutions for eccentricity')
            end
            if nroot >= 2
                n2 = min(n2);
                disp('Error in calculating fluid bearing coefficents - multiple solutions for eccentricity')
            end

            % Compute damping and stiffness matrices
            a = zeros(2,2); b = zeros(2,2);
            n = sqrt(n2);
            eccentricity(ii) = n;
            q1 = (1-n2); q2 = (1+n2); q3 = (1+2*n2); p2 = pi^2;
            de = (p2*q1+16*n2)^1.5;
            a(1,1) = 4*(p2*(2-n2)+16*n2)/de;
            a(2,2) = 4*(p2*q1*q3+32*n2*q2)/(q1*de);
            a(1,2) = pi*(p2*q1^2-16*n^4)/(n*sqrt(q1)*de);
            a(2,1) = -pi*(p2*q1*q3+32*n2*q2)/(n*sqrt(q1)*de);
            b(1,1) = 2*pi*sqrt(q1)*(p2*q3-16*n2)/(n*de);
            b(2,2) = 2*pi*(p2*q1^2+48*n2)/(n*sqrt(q1)*de);
            b(1,2) = -8*(p2*q3-16*n2)/de;
            b(2,1) = b(1,2);
            
            Kb1 = zeros(4,4); Cb1 = zeros(4,4);

            if Linear == 0
                Kb1(1:2,1:2) = (F/c)*a;
                Cb1(1:2,1:2) = (F/(c*rotorSpeed))*b;
            end

        case 8 %Seal 
            if isfield(Rotor.bearing(ii), 'P'), P = Rotor.bearing(ii).P; end  % pressure difference across seal %[Pa]
            if isfield(Rotor.bearing(ii), 'R'), R = Rotor.bearing(ii).R; end % seal radius [m]
            if isfield(Rotor.bearing(ii), 'L'), L = Rotor.bearing(ii).L; end % seal length [m]
            if isfield(Rotor.bearing(ii), 'c'), c = Rotor.bearing(ii).c; end % seal radial clearance [m]
            if isfield(Rotor.bearing(ii), 'V'), V = Rotor.bearing(ii).V; end % seal average axial stream velocity [m/s]
            if isfield(Rotor.bearing(ii), 'fric'), fric = Rotor.bearing(ii).fric; end % friction coefficient
            
            T = L/V;
            sigma = fric*L/c;
            epsilon = pi*sigma*R*P/(6*fric*(1.5+2*sigma));
            mu_0 = 9*sigma/(1.5+2*sigma);
            mu_1 = ( (3+2*sigma)^2*(1.5+2*sigma) - 9*sigma ) / (1.5+2*sigma)^2;
            mu_2 = ( 19*sigma + 18*sigma^2 + 8*sigma^3 ) / (1.5+2*sigma)^3;

            Kb1(1:2,1:2) = epsilon*(mu_0-mu_2*T^2*rotorSpeed^2/4)*eye(2,2) + epsilon*(mu_1*T*rotorSpeed/2)*[0 -1; 1 0];
            Cb1(1:2,1:2) = epsilon*mu_1*T*eye(2,2) + epsilon*(mu_2*T^2*rotorSpeed)*[0 -1; 1 0];
            Mb1(1:2,1:2) = epsilon*mu_2*T^2*eye(2,2);
        case 9      % Ball bearing - simplified model 2DOF
            if isfield(Rotor.bearing(ii), 'n_ball'), n_ball = Rotor.bearing(ii).n_ball; end  % number of balls
            if isfield(Rotor.bearing(ii), 'd_ball'), d_ball = Rotor.bearing(ii).d_ball; end % ball diameter [m]
            if isfield(Rotor.bearing(ii), 'fs'), fs = Rotor.bearing(ii).fs; end % Static Load [N]
            if isfield(Rotor.bearing(ii), 'alpha'), alpha = Rotor.bearing(ii).alpha; end % contact angle [rad]
            
            kb = 13e6;  %[(N^2/m^4)^(1/3)]

            if ~(n_ball == 8 || n_ball == 12 || n_ball == 16)
                error('Bearing Type 9 - Number of balls not valid')
            end
            
            switch n_ball
                case 8
                    coeff = 0.46;
                case 12
                    coeff = 0.64;
                case 16
                    coeff = 0.73;
            end

            kvv = kb * (n_ball^2 * d_ball * fs * cos(alpha)^5)^(1/3);
            kuu = coeff * kvv;
            
            Kb1 = zeros(4,4); Kb1(1,1) = kuu; Kb1(2,2) = kvv;
        case 10 %Squeeze film Dampers with springs
            if isfield(Rotor.bearing(ii), 'eta'), eta = Rotor.bearing(ii).eta; end  % Oil viscosity [Ns/m^2]
            if isfield(Rotor.bearing(ii), 'R'), R = Rotor.bearing(ii).R; end % Damper radius [m]
            if isfield(Rotor.bearing(ii), 'L'), L = Rotor.bearing(ii).L; end % Damper length [m]
            if isfield(Rotor.bearing(ii), 'c'), c = Rotor.bearing(ii).c; end % Damper clearance [m]

            if isfield(Rotor.bearing(ii), 'ks') % Spring stiffness if present [N/m]
                ks = Rotor.bearing(ii).ks;
            else
                ks = 0;
            end 


            Kb1 = zeros(4,4); Cb1 = zeros(4,4);

            Csdf = pi*eta*R*L^3/2/c^3;

            Kb1(1,1) = ks; Kb1(2,2) = ks;
            Cb1(1,1) = Csdf; Cb1(2,2) = Csdf;
        case 11 % Thomas-Alford Forces : Steam & Gas Turbines
            if isfield(Rotor.bearing(ii), 'beta'), beta = Rotor.bearing(ii).beta; end  % Empirical Parameter
            if isfield(Rotor.bearing(ii), 'T'), T = Rotor.bearing(ii).T; end % Turbine Stage Torque [Nm]
            if isfield(Rotor.bearing(ii), 'D'), D = Rotor.bearing(ii).D; end % Blade mean diameter [m]
            if isfield(Rotor.bearing(ii), 'L'), L = Rotor.bearing(ii).L; end % Blade radial length [m]
            
            Kb1 = zeros(4,4); Cb1 = zeros(4,4);

            Ksw = beta*T/D/L;
            Kb1(1,2) = Ksw; Kb1(2,1) = -Ksw;
        case 12 % Scharrer - Childs Labyrinth seals
            if isfield(Rotor.bearing(ii), 'sealType'), sealType = Rotor.bearing(ii).sealType; end
            if isfield(Rotor.bearing(ii), 'P_reserv'), P_reserv = Rotor.bearing(ii).P_reserv; end 
            if isfield(Rotor.bearing(ii), 'P_sump'), P_sump = Rotor.bearing(ii).P_sump; end 
            if isfield(Rotor.bearing(ii), 'Nt'), Nt = Rotor.bearing(ii).Nt; end 
            if isfield(Rotor.bearing(ii), 'Cr'), Cr = Rotor.bearing(ii).Cr; end
            if isfield(Rotor.bearing(ii), 'L'), L = Rotor.bearing(ii).L; end
            if isfield(Rotor.bearing(ii), 'T'), T = Rotor.bearing(ii).T; end
            if isfield(Rotor.bearing(ii), 'Rs'), Rs = Rotor.bearing(ii).Rs; end
            if isfield(Rotor.bearing(ii), 'B'), B = Rotor.bearing(ii).B; end
            if isfield(Rotor.bearing(ii), 'U_inlet'), U_inlet = Rotor.bearing(ii).U_inlet; end
            if isfield(Rotor.bearing(ii), 'omega'), omega = Rotor.bearing(ii).omega; end
            if isfield(Rotor.bearing(ii), 'nu'), nu = Rotor.bearing(ii).nu; end
            
            [~,Kseal,Cseal] = labySeals(sealType,P_reserv,P_sump,...
                Nt,Cr,L,T,Rs,B,U_inlet,omega,nu);

            Kb1([1,2],[1,2]) = Kseal;
            Cb1([1,2],[1,2]) = Cseal;
        case 13 % Stiffness & Damping dependent on RotorSpeed
            Kb1 = zeros(4,4); Cb1 = zeros(4,4);
            % horz vectors only
            if isfield(Rotor.bearing(ii), 'Speed_points'), Speed_points = Rotor.bearing(ii).Speed_points; end
    
            if isfield(Rotor.bearing(ii), 'kxx_points'), kxx_points = Rotor.bearing(ii).kxx_points; end
            if isfield(Rotor.bearing(ii), 'kxy_points'), kxy_points = Rotor.bearing(ii).kxy_points; end
            if isfield(Rotor.bearing(ii), 'kyx_points'), kyx_points = Rotor.bearing(ii).kyx_points; end
            if isfield(Rotor.bearing(ii), 'kyy_points'), kyy_points = Rotor.bearing(ii).kyy_points; end

            if isfield(Rotor.bearing(ii), 'cxx_points'), cxx_points = Rotor.bearing(ii).cxx_points; end
            if isfield(Rotor.bearing(ii), 'cxy_points'), cxy_points = Rotor.bearing(ii).cxy_points; end
            if isfield(Rotor.bearing(ii), 'cyx_points'), cyx_points = Rotor.bearing(ii).cyx_points; end
            if isfield(Rotor.bearing(ii), 'cyy_points'), cyy_points = Rotor.bearing(ii).cyy_points; end
            
            K_table = [kxx_points(:), kxy_points(:), kyx_points(:), kyy_points(:)];
            C_table = [cxx_points(:), cxy_points(:), cyx_points(:), cyy_points(:)];

            K_vec = interp1(Speed_points, K_table, rotorSpeed, 'linear', 'extrap');
            C_vec = interp1(Speed_points, C_table, rotorSpeed, 'linear', 'extrap');
            
            Kb1([1,2],[1,2]) = [K_vec(1), K_vec(2); K_vec(3), K_vec(4)];
            Cb1([1,2],[1,2]) = [C_vec(1), C_vec(2); C_vec(3), C_vec(4)];
    end
    
    nnode = bearings(ii).node;
    dof = (4*nnode-3):4*nnode;
    Kb(dof,dof) = Kb(dof,dof) + Kb1;
    Cb(dof,dof) = Cb(dof,dof) + Cb1;
    Mb(dof,dof) = Mb(dof,dof) + Mb1;
end
end

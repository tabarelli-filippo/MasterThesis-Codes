function [critical_speeds,mode_shapes] = critSpeeds(Rotor,NX,isDamped,method,options)
% CRITSPEEDS  Computes critical speeds of a rotor-bearing system.
%
%   A critical speed is the shaft speed at which the excitation frequency
%   (NX * rotorSpeed) coincides with a natural frequency of the system.
%   Three solution methods are available:
%
%   Method 1 - Direct: solves a single eigenvalue problem at the target
%     speed, treating the system matrices as evaluated at Omega such that
%     the NX excitation is synchronous. Efficient for simple cases.
%
%   Method 2 - Iterative (no initial guess): starts from a nominal speed
%     estimate and iterates each critical speed separately by updating the
%     speed-dependent matrices until convergence. No user guess needed.
%
%   Method 3 - Iterative (with initial guesses): same as Method 2 but
%     uses user-supplied initial estimates to target specific critical
%     speeds. Recommended when critical speeds are close together.
%
% SYNTAX
%   critical_speeds = critSpeeds(Rotor, NX, isDamped, method)
%   [critical_speeds, mode_shapes] = critSpeeds(Rotor, NX, isDamped, method, ...
%       "num_crit", nc, "max_iter", mi, "toll", tol, ...
%       "initial_estimates", init_est)
%
% INPUT ARGUMENTS
%   Rotor       - (1x1 struct) Rotor data structure (see rotorMatrix,
%                 bearingMatrix for field definitions)
%   NX          - (double, scalar) Engine order of the excitation, i.e.
%                 the ratio omega_excitation / rotorSpeed. For synchronous
%                 excitation (mass unbalance) use NX = 1.
%   isDamped    - (logical) false → use undamped natural frequencies;
%                            true  → use damped natural frequencies
%   method      - (double) Solution method selector: 1, 2, or 3 (see above)
%
% NAME-VALUE OPTIONS (methods 2 and 3)
%   "num_crit"          - (positive integer, default 5) Number of critical
%                         speeds to compute
%   "max_iter"          - (positive integer, default 20) Maximum number of
%                         iterations per critical speed
%   "toll"              - (positive double, default 1e-3) Relative
%                         convergence tolerance: |Omega_new - Omega_old|
%                         / Omega_old < toll
%   "initial_estimates" - (vector double) Initial speed estimates [rad/s]
%                         for each critical speed (method 3 only)
%
% OUTPUT ARGUMENTS
%   critical_speeds - (num_crit x 1 double) Critical speeds [rad/s].
%                     For method 1, all critical speeds up to num_crit are
%                     returned simultaneously.
%   mode_shapes     - (ndof x num_crit double, optional) Mode shape vectors
%                     at each critical speed. Returned only when two output
%                     arguments are requested.
%
% NOTES
%   - Fluid-film bearings (type 7) require rotorSpeed > 0; a dummy speed
%     of 100 rad/s is used during matrix initialisation.
%   - A warning is issued if an iterative method fails to converge within
%     max_iter iterations.
%   - ndof = 4 * n_nodes (global degrees of freedom)
%
% EXAMPLE
%   % Synchronous critical speeds, iterative method
%   crit = critSpeeds(Rotor, 1, false, 2, "num_crit", 3);
%   fprintf('1st critical: %.1f rpm\n', crit(1)*60/(2*pi));
%
%   % With initial guesses
%   init = [300; 800; 1500];   % [rad/s]
%   [crit, shapes] = critSpeeds(Rotor, 1, true, 3, ...
%       "num_crit", 3, "initial_estimates", init);
%
% SEE ALSO
%   charRoots, rotorMatrix, bearingMatrix, FRF

arguments (Input)
    Rotor (1,1) struct
    NX (1,1) double
    isDamped (1,1) logical
    method (1,1) double

    options.num_crit(1,1) double {mustBeInteger, mustBePositive} = 5
    options.max_iter(1,1) double {mustBeInteger, mustBePositive} = 20
    options.toll (1,1) double {mustBePositive} = 1e-3
    options.initial_estimates (:,:) double {mustBeVector}
end

rotorSpeed = 100; % dummy speed, some bearings are not defined at 0 [rad/s]
n_nodes = numel(Rotor.nodes);
ndof = 4 * n_nodes;

% define global matrices
[M0,C0,C1,K0,K1] = rotorMatrix(Rotor);
[Mb,Cb,Kb,zero_dof,~] = bearingMatrix(Rotor,rotorSpeed);

% sort out zeroed DoF
dof = 1:ndof;
dof(zero_dof) = [];
nzero = length(zero_dof);
ncdof = ndof - nzero;

eigenvalues = zeros(2*ncdof,1);
%% direct method
if method == 1
    % compute global matrices
    j = 1i;
    MM = -(NX^2)*(M0 + Mb) + j*NX*C1;
    CC = NX*(C0 + Cb) -j*K1;
    KK = K0 + Kb;
    AA = [zeros(ncdof,ncdof) eye(ncdof,ncdof); -MM(dof,dof)\KK(dof,dof) -MM(dof,dof)\CC(dof,dof)];
    
    if nargout<=1   %eigs only
        eigenvalues = sort(eig(AA));
    end

    if nargout==2   %eigs and eigvectors
        [eigenvectors,eigenvalues] = eig(AA);
        [eigenvalues,isort] = sort(diag(eigenvalues));
        eigenvectors = eigenvectors(:,isort);
        mode_shapes = zeros(ndof,ncdof);
        mode_shapes(dof,:) = eigenvectors(1:ncdof,(1:2:2*ncdof));       
    end
    
    if isDamped % natural damped frequencies
        critical_speeds = abs(real(eigenvalues(1:2:2*ncdof)));
    else % natural frequencies
        critical_speeds = abs(eigenvalues(1:2:2*ncdof));
    end
end 

%% iterative method - no guesses
if method == 2
    number_criticals = options.num_crit;
    max_iterations = options.max_iter;
    convergence_tol = options.toll;

    % Initial guess
    rotorSpeed_initial = pi*50/3; %[rad/s]
    
    % compute global matrices - speed dependent
    M = M0 + Mb;
    K = K0 + Kb + rotorSpeed_initial*K1;
    C = C0 + Cb + rotorSpeed_initial*C1;
    % initial eigenvalues
    AA = [zeros(ncdof,ncdof) eye(ncdof,ncdof); -M(dof,dof)\K(dof,dof) -M(dof,dof)\C(dof,dof)];
    eigenvalues_initial = sort(eig(AA));

    critical_speeds = zeros(number_criticals,1);
    for icritical = 1:number_criticals
        % initialization
        ieig = 2*icritical-1;
        rel_change = 1;
        iteration = 0;
        
        if isDamped
            critical_i = abs(imag(eigenvalues_initial(ieig)))/NX;
        else
            critical_i = abs(eigenvalues_initial(ieig))/NX;
        end

        while (rel_change>convergence_tol) && (iteration<max_iterations)
            %define new matrices
            [Mb,Cb,Kb,~,~] = bearingMatrix(Rotor,critical_i);
            M = M0 + Mb;
            K = K0 + Kb + critical_i*K1;
            C = C0 + Cb + critical_i*C1;
            %new critical speeds
            AA = [zeros(ncdof,ncdof) eye(ncdof,ncdof); -M(dof,dof)\K(dof,dof) -M(dof,dof)\C(dof,dof)];
            eigenvalues = sort(eig(AA));
            critical_i_old = critical_i;

            if isDamped
                critical_i = abs(imag(eigenvalues(ieig)))/NX;
            else
                critical_i = abs(eigenvalues(ieig))/NX;
            end

            if critical_i_old == 0
                rel_change = 2*convergence_tol;
            else
                rel_change = abs( (critical_i - critical_i_old)/critical_i_old );
            end
            iteration = iteration + 1;
        end

        if (rel_change > convergence_tol)
        warning('CRITSPEEDS:NoConvergence', 'Convergence for Critical Speed %d not reached after %d iterations', ...
                icritical, (iteration-1));
        end
        critical_speeds(icritical) = critical_i;
    end
    
    % if eigenvectors are requested
    if nargout == 2
        mode_shapes = zeros(ndof,number_criticals);
        for icritical = 1:number_criticals
            [Mb,Cb,Kb,~,~] = bearingMatrix(Rotor,critical_speeds(icritical));
            M = M0 + Mb;
            K = K0 + Kb + critical_speeds(icritical)*K1;
            C = C0 + Cb + critical_speeds(icritical)*C1;
            AA = [zeros(ncdof,ncdof) eye(ncdof,ncdof); -M(dof,dof)\K(dof,dof) -M(dof,dof)\C(dof,dof)];
            [eigenvectors_i,eigenvalues] = eig(AA);
            [~,isort] = sort(diag(eigenvalues));
            eigenvectors_i = eigenvectors_i(:,isort);
            mode_shapes(dof,icritical) = eigenvectors_i(1:ncdof,2*icritical-1);
        end
    end
end

%% iterative method - initial guesses
if method == 3
    number_criticals = options.num_crit;
    max_iterations = options.max_iter;
    convergence_tol = options.toll;
    initial_estimates = options.initial_estimates;

    % Initial guess
    rotorSpeed_initial = pi*50/3; %[rad/s]
    
    critical_speeds = zeros(number_criticals,1);
    for icritical = 1:number_criticals
        rel_change = 1;
        iteration = 0;
        critical_i = abs(initial_estimates(icritical));

        while (rel_change>convergence_tol) && (iteration<max_iterations)
            %define new matrices
            [Mb,Cb,Kb,~,~] = bearingMatrix(Rotor,critical_i);
            M = M0 + Mb;
            K = K0 + Kb + critical_i*K1;
            C = C0 + Cb + critical_i*C1;
            %new critical speeds
            AA = [zeros(ncdof,ncdof) eye(ncdof,ncdof); -M(dof,dof)\K(dof,dof) -M(dof,dof)\C(dof,dof)];
            eigenvalues = sort(eig(AA));
            critical_i_old = critical_i;

            if isDamped
                critical_est = abs(imag(eigenvalues))/NX;
            else
                critical_est = abs(eigenvalues)/NX;
            end

            [~,index_c] = min( abs( critical_est - initial_estimates(icritical) ) );
            critical_i = critical_est(index_c);
            
            if critical_i_old == 0
                rel_change = 2*convergence_tol;
            else
                rel_change = abs( (critical_i - critical_i_old)/critical_i_old );
            end

            iteration = iteration + 1;
        end

        if (rel_change > convergence_tol)
        warning('CRITSPEEDS:NoConvergence', 'Convergence for Critical Speed %d not reached after %d iterations', ...
                icritical, (iteration-1));
        end
        critical_speeds(icritical) = critical_i;
    end
    
    % if eigenvectors are requested
    if nargout == 2
        mode_shapes = zeros(ndof,number_criticals);
        for icritical = 1:number_criticals
            [Mb,Cb,Kb,~] = bearingMatrix(Rotor,critical_speeds(icritical));
            M = M0 + Mb;
            K = K0 + Kb + critical_speeds(icritical)*K1;
            C = C0 + Cb + critical_speeds(icritical)*C1;
            AA = [zeros(ncdof,ncdof) eye(ncdof,ncdof); -M(dof,dof)\K(dof,dof) -M(dof,dof)\C(dof,dof)];
            [eigenvectors_i,eigenvalues] = eig(AA);
            [eigenvalues,isort] = sort(diag(eigenvalues));
            eigenvectors_i = eigenvectors_i(:,isort);

            if isDamped
                critical_est = abs(imag(eigenvalues))/NX;
            else
                critical_est = abs(eigenvalues)/NX;
            end

            [~,index_c] = min( abs( critical_est - critical_speeds(icritical) ) );           
            mode_shapes(dof,icritical) = eigenvectors_i(1:ncdof,index_c);
        end
    end
end
end

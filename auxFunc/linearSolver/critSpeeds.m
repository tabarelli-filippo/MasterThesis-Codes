function [critical_speeds,mode_shapes] = critSpeeds(Rotor,NX,isDamped,method,options)
%CRITSPEEDS evaluates critical speeds for a rotor, using direct approach or
%   iterative method.
%INPUT: Rotor       Structure - Rotor
%
%       NX          double
%
%       isDamped    parameter to plot undamped natural frequencies (default,
%                   dampedFreq=0) or undamped natural frequencies (dampedFreq=1)
%       
%       method      method=1 - Direct Method
%                   method=2 - Iterative Method
%
%   [critical_speeds, mode_shapes] = critSpeeds(Rotor, NX, isDamped, method, ...
%                              "num_crit", num_crit,"max_iter", max_iter, "toll", toll,...
%                              "initial_estimates", init_est);
%
%   num_crit = amount of critical speeds needed
%   max_iter = maximum interations number
%   toll = convergence tollerance
%   initial_estimates

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
    % initial eigevalues
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
    
    % if eigevectors are requested
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
    
    % if eigevectors are requested
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
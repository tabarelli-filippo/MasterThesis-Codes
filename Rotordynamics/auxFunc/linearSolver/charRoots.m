function [eigenvalues,eigenvectors,kappa,mass_Matrix] = charRoots(Rotor,rotorSpeed)
% CHARROOTS  Computes eigenvalues, eigenvectors, and whirl ratios of a
%   rotor-bearing system over a range of rotational speeds.
%
%   The function builds the speed-dependent state-space matrix of the rotor
%   and solves the associated eigenvalue problem at each requested speed.
%   Overdamped and spurious (non-oscillatory) modes are automatically
%   filtered out. Modes are sorted by ascending damped natural frequency.
%
%   The state-space formulation is:
%       A(Omega) = [  0        I   ]
%                  [ -M\K   -M\C  ]
%   where M, K, C include contributions from shaft, disks, and bearings,
%   and the gyroscopic term is speed-proportional.
%
% SYNTAX
%   eigenvalues                        = charRoots(Rotor, rotorSpeed)
%   [eigenvalues, eigenvectors]        = charRoots(Rotor, rotorSpeed)
%   [eigenvalues, eigenvectors, kappa] = charRoots(Rotor, rotorSpeed)
%   [eigenvalues, eigenvectors, kappa, mass_Matrix] = charRoots(Rotor, rotorSpeed)
%
% INPUT ARGUMENTS
%   Rotor       - (1x1 struct) Rotor data structure containing node, shaft,
%                 disk, bearing, and forcing definitions (see rotorMatrix,
%                 bearingMatrix)
%   rotorSpeed  - (vector, double) Array of rotational speeds [rad/s] at
%                 which eigenvalues are computed
%
% OUTPUT ARGUMENTS
%   eigenvalues  - (ncdof x nspeed complex double) Eigenvalues of the
%                  state-space system at each speed. Each eigenvalue is of
%                  the form sigma + i*omega_d, where sigma is the growth/
%                  decay rate [1/s] and omega_d is the damped natural
%                  frequency [rad/s]. Entries may be NaN if fewer valid
%                  modes than ncdof are found at a given speed.
%   eigenvectors - (ndof x ncdof [x nspeed] complex double) Physical part
%                  of the state-space eigenvectors (displacement partition),
%                  zero-padded at constrained DOFs. 3-D array if nspeed > 1.
%   kappa        - (ndof x ncdof [x nspeed] double) Whirl ratio at each
%                  node for each mode and speed, computed via modeWhirl.
%                  kappa > 0 indicates forward whirl (FW);
%                  kappa < 0 indicates backward whirl (BW).
%                  3-D array if nspeed > 1.
%   mass_Matrix  - (ndof x ndof double) Global mass matrix M0 (shaft and
%                  disks only, bearing mass not included)
%
% NOTES
%   - Only eigenvalues with Im(lambda) < -1e-1 rad/s are retained
%     (oscillatory modes with negative imaginary part per the sign
%     convention used in this solver).
%   - Constrained DOFs (from bearingMatrix) are removed before solving.
%   - ndof  = 4 * n_nodes  (global degrees of freedom)
%   - ncdof = ndof - nzero (constrained DOFs removed)
%
% EXAMPLE
%   speeds = linspace(100, 1000, 50);   % [rad/s]
%   [evals, evecs, kap] = charRoots(Rotor, speeds);
%   % Campbell diagram: plot imaginary parts vs speeds
%   figure; plot(speeds, abs(imag(evals)) / (2*pi))
%   xlabel('Speed [rad/s]'); ylabel('Frequency [Hz]')
%
% SEE ALSO
%   rotorMatrix, bearingMatrix, modeWhirl, sortModesMAC, critSpeeds

arguments (Input)
    Rotor (1,1) struct
    rotorSpeed (:,:) double {mustBeVector}
end

if nargout > 1
    % initialize solutions
    n_nodes = numel(Rotor.nodes);
    nspeed = length(rotorSpeed);
    
    ndof = 4 * n_nodes;
    % define global matrices
    [M0,C0,C1,K0,K1] = rotorMatrix(Rotor);
    [Mb,Cb,Kb,zero_dof,~] = bearingMatrix(Rotor,rotorSpeed(1));
    
    % sort out zeroed DoF
    dof = 1:ndof;
    dof(zero_dof) = [];
    nzero = length(zero_dof);
    ncdof = ndof - nzero;
    
    % initialize to NaN to handle filtered modes
    eigenvalues = nan(ncdof,nspeed);
    
    if nspeed==1
        eigenvectors = zeros(ndof,ncdof);
    else
        eigenvectors = zeros(ndof,ncdof,nspeed);
    end
    
    % finding eigenvalues for each speed
    for ii = 1:nspeed
        [Mb,Cb,Kb,~,~] = bearingMatrix(Rotor,rotorSpeed(ii));
        % speed-dependent global matrices
        M = M0 + Mb;
        K = K0 + Kb + rotorSpeed(ii)*K1;
        C = C0 + Cb + rotorSpeed(ii)*C1;
        
        % State Space Matrix
        AA = [zeros(ncdof,ncdof) eye(ncdof,ncdof); -M(dof,dof)\K(dof,dof) -M(dof,dof)\C(dof,dof)];
        
        [V_raw, D_raw] = eig(AA);
        vals_raw = diag(D_raw);

        vibration_threshold = 1e-1; % Rad/s
        
        idx_valid = find(imag(vals_raw) < -vibration_threshold);
        
        vals_filtered = vals_raw(idx_valid);
        vecs_filtered = V_raw(:, idx_valid);
        
        [~, sort_idx] = sort(abs(imag(vals_filtered)));
        
        vals_sorted = vals_filtered(sort_idx);
        vecs_sorted = vecs_filtered(:, sort_idx);
        
        n_found = min(length(vals_sorted), ncdof);
        
        eigenvalues(1:n_found, ii) = vals_sorted(1:n_found);
        physical_modes = vecs_sorted(1:ncdof, 1:n_found);
        
        if nspeed == 1
            eigenvectors(dof, 1:n_found) = physical_modes; 
        else
            eigenvectors(dof, 1:n_found, ii) = physical_modes; 
        end
    end
end

if nargout > 2
    if nspeed==1
        kappa = zeros(ndof,ncdof);
    else
        kappa = zeros(ndof,ncdof,nspeed);
    end
    
    for i = 1:nspeed        % for each speed
        valid_modes_count = sum(~isnan(eigenvalues(:,i)));
        
        for id = 1:valid_modes_count  % for each valid eigenvector
            if nspeed == 1
                evector = eigenvectors(:,id);
            else
                evector = eigenvectors(:,id,i);
            end
            
            kk = modeWhirl(evector(1:2:ndof),evector(2:2:ndof));
            
            if imag(eigenvalues(id,i)) < 0
                kk = -kk;
            end
            
            if nspeed == 1
                kappa(1:2:ndof,id) = kk;
                kappa(2:2:ndof,id) = kk;
            else
                kappa(1:2:ndof,id,i) = kk;
                kappa(2:2:ndof,id,i) = kk;
            end
        end
    end
end

if nargout > 3
    mass_Matrix = M0;
end

end

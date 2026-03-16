function [eigenvalues,eigenvectors,kappa,mass_Matrix] = charRoots(Rotor,rotorSpeed)
%CHARROOTS: computes Natural frequencies and Modes of a rotor.
% Includes filtering for overdamped/spurious modes.
%INPUT: Rotor       Structure of the rotor
%       rotorSpeed  vector of speeds to evaluate

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
    
    % MODIFICA 1: Inizializziamo a NaN per gestire i modi filtrati
    eigenvalues = nan(ncdof,nspeed);
    
    if nspeed==1
        eigenvectors = zeros(ndof,ncdof);
    else
        eigenvectors = zeros(ndof,ncdof,nspeed);
    end
    
    % finding eigevalues for each speed
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
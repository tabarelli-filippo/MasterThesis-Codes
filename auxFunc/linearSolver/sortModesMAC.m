function [sorted_vals, sorted_vecs, sorted_kappa] = sortModesMAC(eigs, eigenvectors, kappa)
% SORTMODESMAC modes sorting using MAC coefficient
% INPUTS:   eigs   : eigenvalues matric
%           eigenvectors: eigenvector matrix
%           kappa : whirl matrix

[~, n_modes, n_steps] = size(eigenvectors);
%% initialize
sorted_vals = eigs;
sorted_vecs = eigenvectors;

has_kappa = (nargin > 2 && ~isempty(kappa));
if has_kappa
    sorted_kappa = kappa;
else
    sorted_kappa = [];
end

%% reference velocity: Omega = 0;
ref_phi = eigenvectors(:, :, 1);
for k = 2:n_steps
    curr_phi = eigenvectors(:, :, k);
    curr_val = eigs(:, k);

    if has_kappa
        if ndims(kappa) == 3
            curr_kap = kappa(:, :, k); % Tensore
        else
            curr_kap = kappa(:, k);    % Matrice
        end
    end

    %% MAC evaluation
    cross_prod = abs(ref_phi' * curr_phi).^2;
    auto_ref   = diag(ref_phi' * ref_phi);
    auto_curr  = diag(curr_phi' * curr_phi);

    denom = auto_ref * auto_curr';
    MAC_matrix = cross_prod ./ denom;

    new_order = zeros(1, n_modes);
    temp_MAC = MAC_matrix;

    for i = 1:n_modes
        % search for matching
        [~, best_match_idx] = max(temp_MAC(i, :));
        
        new_order(i) = best_match_idx;
        % eliminating already used column
        temp_MAC(:, best_match_idx) = -1;
    end

    sorted_vals(:, k)     = curr_val(new_order);
    sorted_vecs(:, :, k)  = curr_phi(:, new_order);

    if has_kappa
        if ndims(kappa) == 3
            sorted_kappa(:, :, k) = curr_kap(:, new_order);
        else
            sorted_kappa(:, k) = curr_kap(new_order); 
        end
    end
    ref_phi = sorted_vecs(:, :, k);
end
end
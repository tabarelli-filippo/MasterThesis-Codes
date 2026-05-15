function [sorted_vals, sorted_vecs, sorted_kappa] = sortModesMAC(eigs, eigenvectors, kappa)
% SORTMODESMAC  Sorts rotor modes across speed steps using the Modal
%   Assurance Criterion (MAC) to ensure mode tracking consistency along
%   a Campbell diagram or speed sweep.
%
%   At each speed step, each mode of the current step is matched to the
%   mode of the previous step that maximises the MAC value. This prevents
%   mode-crossing artefacts that would otherwise appear when eigenvalues
%   approach each other or swap order.
%
%   The MAC between two complex vectors phi_i and phi_j is defined as:
%       MAC(i,j) = |phi_i' * phi_j|^2 / (|phi_i'*phi_i| * |phi_j'*phi_j|)
%
%   A greedy assignment strategy is used: for each reference mode (row),
%   the current mode with the highest MAC is selected, and that column is
%   removed from consideration for subsequent assignments.
%
% SYNTAX
%   [sorted_vals, sorted_vecs]               = sortModesMAC(eigs, eigenvectors)
%   [sorted_vals, sorted_vecs, sorted_kappa] = sortModesMAC(eigs, eigenvectors, kappa)
%
% INPUT ARGUMENTS
%   eigs         - (ncdof x nspeed complex double) Eigenvalue matrix as
%                  returned by charRoots
%   eigenvectors - (ndof x ncdof x nspeed complex double) Eigenvector
%                  tensor as returned by charRoots
%   kappa        - (ndof x ncdof x nspeed double, optional) Whirl ratio
%                  tensor as returned by charRoots. If omitted or empty,
%                  sorted_kappa is returned as [].
%
% OUTPUT ARGUMENTS
%   sorted_vals  - (ncdof x nspeed complex double) Eigenvalues reordered
%                  by MAC-based mode tracking at each speed step
%   sorted_vecs  - (ndof x ncdof x nspeed complex double) Corresponding
%                  reordered eigenvectors
%   sorted_kappa - (ndof x ncdof x nspeed double) Corresponding reordered
%                  whirl ratios (empty if kappa was not supplied)
%
% NOTES
%   - The first speed step (index 1) is always used as the reference and
%     is returned unchanged.
%   - The reference is updated after each step (sliding-window tracking),
%     which improves robustness for large speed ranges.
%   - The greedy MAC assignment may fail when two modes have nearly
%     identical shapes; in such cases, manual inspection of the Campbell
%     diagram is recommended.
%
% EXAMPLE
%   speeds = linspace(0, 1500, 100);
%   [evals, evecs, kap] = charRoots(Rotor, speeds);
%   [evals_s, evecs_s, kap_s] = sortModesMAC(evals, evecs, kap);
%   % evals_s now contains continuously tracked modal branches
%
% SEE ALSO
%   charRoots, modeWhirl

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
            curr_kap = kappa(:, :, k); % 3-D tensor
        else
            curr_kap = kappa(:, k);    % 2-D matrix
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
        % search for best matching mode (greedy assignment)
        [~, best_match_idx] = max(temp_MAC(i, :));
        
        new_order(i) = best_match_idx;
        % remove already assigned column to avoid duplicate matching
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

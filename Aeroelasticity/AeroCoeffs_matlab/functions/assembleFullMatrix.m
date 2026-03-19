function [M_full, K_full] = assembleFullMatrix(M1, K1, M2, K2, N)
%ASSEMBLE_FULL_ROTOR full matrix from sector
%   INPUT:
%   - M1, K1:   sector matrices  [Dofs x Dofs]
%   - M2, K2:   coupling matrices [Dofs x Dofs]
%   - N: number of blades
%
%   OUTPUT:
%   - M_full, K_full: Sparse full global matrices [N*Dofs x N*Dofs]

    K1  = sparse(K1);
    K2 = sparse(K2);
    M1  = sparse(M1);
    M2 = sparse(M2);

    %% 1. Circular matrices permutations
    I_N = sparse(eye(N));
    
    % P_prev: Lower Diag + upper right corner (Wrap-around)
    % [ 0 0 ... 1 ]x
    % [ 1 0 ... 0 ]
    % [ 0 1 ... 0 ]
    P_prev = sparse(diag(ones(N-1,1), -1));
    P_prev(1, end) = 1; 
    
    % P_next: Upper Diag + lower left corner (Wrap-around)
    P_next = P_prev';

    %% 2. Assembling
    K_full = kron(I_N, K1) + kron(P_prev, K2) + kron(P_next, K2');

    %% 3. Assemblaggio Matrice Massa
    M_full = kron(I_N, M1) +kron(P_prev, M2) + kron(P_next, M2');
end
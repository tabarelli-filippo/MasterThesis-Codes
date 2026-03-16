function delta = LogDecrement(eigenvalues)
%LOGDECREMENT: evaluates log decrement for eigenvalues
%
%INPUT: eigenvalues
    arguments
        eigenvalues double
    end

    sigma = real(eigenvalues);
    omega_d = abs(imag(eigenvalues)); 
    omega_threshold = 1e-10;
    
    delta = nan(size(eigenvalues));
    
    is_oscillatory = omega_d > omega_threshold;
    
    delta(is_oscillatory) = -2 * pi * (sigma(is_oscillatory) ./ omega_d(is_oscillatory));
end
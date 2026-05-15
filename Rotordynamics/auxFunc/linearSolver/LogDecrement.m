function delta = LogDecrement(eigenvalues)
% LOGDECREMENT  Computes the logarithmic decrement for a set of complex
%   eigenvalues from a rotor-bearing system.
%
%   The logarithmic decrement is a dimensionless measure of damping defined
%   as the natural logarithm of the ratio of successive amplitudes in a
%   free oscillation. For a complex eigenvalue lambda = sigma + i*omega_d:
%
%       delta = -2 * pi * sigma / omega_d
%
%   where sigma = Re(lambda) is the decay rate [1/s] and
%   omega_d = |Im(lambda)| is the damped natural frequency [rad/s].
%
%   A positive delta indicates stable (decaying) oscillation.
%   delta = 0 corresponds to undamped vibration.
%   Negative delta indicates an unstable (growing) oscillation.
%
% SYNTAX
%   delta = LogDecrement(eigenvalues)
%
% INPUT ARGUMENTS
%   eigenvalues - (array, complex double) Complex eigenvalues of the
%                 rotor system, as returned by charRoots or eig.
%                 Can be a vector or matrix of arbitrary size.
%
% OUTPUT ARGUMENTS
%   delta - (same size as eigenvalues, double) Logarithmic decrement for
%           each eigenvalue. Entries are NaN for purely real (non-
%           oscillatory, overdamped) eigenvalues, identified by
%           |Im(lambda)| <= 1e-10 rad/s.
%
% NOTES
%   - Non-oscillatory modes (|Im(lambda)| ≤ 1e-10) yield NaN, as the
%     logarithmic decrement is not meaningful for overdamped modes.
%   - The sign convention follows the standard rotor dynamics convention
%     where eigenvalues have negative real parts for stable systems.
%
% EXAMPLE
%   speeds = linspace(100, 2000, 50);
%   [evals, ~] = charRoots(Rotor, speeds);
%   delta = LogDecrement(evals);
%   % Plot logarithmic decrement for the first mode vs speed
%   figure; plot(speeds * 60/(2*pi), delta(1, :));
%   xlabel('Speed [rpm]'); ylabel('Log. decrement [-]');
%   yline(0, '--r'); title('Stability map - Mode 1');
%
% SEE ALSO
%   charRoots, critSpeeds

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

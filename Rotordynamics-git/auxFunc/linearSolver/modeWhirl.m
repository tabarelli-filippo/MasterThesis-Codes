function kappa = modeWhirl(u,v)
% MODEWHIRL  Computes the whirl ratio kappa for each node of a rotor mode
%   shape to determine the direction of precession (forward or backward).
%
%   For each node, the lateral displacement orbit is an ellipse described
%   by the complex amplitudes u (x-direction) and v (y-direction). The
%   whirl ratio kappa is defined as the ratio of the minor to the major
%   semi-axis of this ellipse, with a sign that encodes the whirl direction:
%
%       kappa > 0  →  Forward Whirl (FW): precession in the same direction
%                     as shaft rotation
%       kappa < 0  →  Backward Whirl (BW): precession opposite to rotation
%       |kappa| = 1 →  Circular orbit (purely forward or backward)
%       kappa = 0   →  Planar (degenerate) orbit; direction indeterminate
%
%   The algorithm builds the 2x2 orbit tensor T from the amplitude and
%   phase of u and v, computes H = T*T', and extracts the semi-axis lengths
%   as the square roots of the eigenvalues of H. The sign is determined
%   from the angular difference of the phases of v and u.
%
% SYNTAX
%   kappa = modeWhirl(u, v)
%
% INPUT ARGUMENTS
%   u - (1 x npts complex double) Complex x-displacement amplitude at each
%       node of the mode shape (e.g., extracted from eigenvectors at odd
%       indices: eigenvector(1:2:ndof))
%   v - (1 x npts complex double) Complex y-displacement amplitude at each
%       node of the mode shape (e.g., extracted from eigenvectors at even
%       indices: eigenvector(2:2:ndof))
%
% OUTPUT ARGUMENTS
%   kappa - (1 x npts double) Whirl ratio at each node. Values are in the
%           range [-1, 1]. kappa = 0 at nodes with negligible displacement.
%
% NOTES
%   - Nodes with |u|*|v| < 1e-16 are assigned kappa = 0 to avoid division
%     by zero.
%   - The sign convention assumes that eigenvectors follow the standard
%     complex notation used in charRoots (negative imaginary part for FW
%     modes).
%   - This function is called internally by charRoots to populate the
%     kappa output array.
%
% EXAMPLE
%   [~, evecs, ~] = charRoots(Rotor, 500);
%   mode1 = evecs(:, 1);          % first mode at 500 rad/s
%   u = mode1(1:2:end);           % x-components
%   v = mode1(2:2:end);           % y-components
%   kap = modeWhirl(u, v);
%   % kap(i) > 0: node i precesses forward
%
% SEE ALSO
%   charRoots, sortModesMAC

npts = length(u);
ru = abs(u);
rv = abs(v);
eta_u = angle(u);
eta_v = angle(v);
kappa = zeros(1,npts);

for ii = 1:npts
     if ru(ii)* rv(ii) < 1e-16
        kappa(ii) = 0;
    else
        T = [ru(ii)*cos(eta_u(ii)), -ru(ii)*sin(eta_u(ii)); rv(ii)*cos(eta_v(ii)), -rv(ii)*sin(eta_v(ii))];

        H = T * T';
        L = eigs(H);
        lambda = sqrt(L);
        kappa(ii) = min(lambda) / max(lambda);

        % determine whirl direction from phase difference
        diff = mod(eta_v(ii) - eta_u(ii),2*pi);
        if diff > 0 && diff < pi
            kappa(ii) = - kappa(ii);
        end
     end
end

end

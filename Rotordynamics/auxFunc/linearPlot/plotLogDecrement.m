function [] = plotLogDecrement(rotorSpeed, eigenvalues, options)
% PLOTLOGDECREMENT  Plots the logarithmic decrement versus rotor speed to
%   assess modal stability across the operating range.
%
%   For each eigenvalue lambda = sigma + i*omega_d, the logarithmic
%   decrement is computed as:
%       delta = -2 * pi * sigma / omega_d
%
%   Stable modes (Re(lambda) <= threshold) are plotted as blue circles;
%   unstable modes (Re(lambda) > threshold) are plotted as filled red
%   circles. A horizontal dashed line at delta = 0 marks the stability
%   boundary. Modes with |Im(lambda)| <= 1e-10 (overdamped) are excluded.
%
% SYNTAX
%   plotLogDecrement(rotorSpeed, eigenvalues)
%   plotLogDecrement(rotorSpeed, eigenvalues, "threshold_stab", thr)
%
% INPUT ARGUMENTS
%   rotorSpeed  - (vector, double) Rotor speeds at which eigenvalues were
%                 evaluated [rpm]. Must be consistent with the columns (or
%                 rows) of the eigenvalues matrix.
%   eigenvalues - (matrix or vector, complex double) Complex eigenvalues
%                 as returned by charRoots. Each entry is of the form
%                 sigma + i*omega_d [rad/s]. The function accepts both
%                 (ncdof x nspeed) and (nspeed x ncdof) layouts.
%
% NAME-VALUE OPTIONS
%   "threshold_stab" - (double, default 1e-6) Stability threshold on
%                      Re(lambda) [1/s]. Eigenvalues with Re > threshold
%                      are classified as unstable and highlighted in red.
%
% OUTPUT
%   A figure titled 'Damping vs Rotor Speed' is created with:
%     - Blue circles : stable modes (Re(lambda) <= threshold)
%     - Filled red circles : unstable modes (Re(lambda) > threshold)
%     - Black dashed line : delta = 0 stability boundary
%   The y-axis is limited to the physically meaningful range of delta,
%   clipping outliers with |delta| > 10 for readability.
%
% EXAMPLE
%   speeds_rpm = linspace(500, 5000, 100);
%   speeds_rads = speeds_rpm * pi/30;
%   [evals, ~] = charRoots(Rotor, speeds_rads);
%   plotLogDecrement(speeds_rpm, evals);
%
% SEE ALSO
%   LogDecrement, charRoots, plotRootLocus, plotCampbell

arguments (Input)
    rotorSpeed (:,:) double {mustBeVector}
    eigenvalues (:,:) double
    options.threshold_stab (1,1) double = 1e-6
end

threshold = options.threshold_stab;


re_eigs = real(eigenvalues);
im_eigs = abs(imag(eigenvalues));

is_oscillatory = im_eigs > 1e-10; 
log_dec = zeros(size(eigenvalues)); 
log_dec(is_oscillatory) = -2 * pi * re_eigs(is_oscillatory) ./ im_eigs(is_oscillatory);
log_dec(~is_oscillatory) = NaN; 

stable_idx = re_eigs <= threshold;
unstable_idx = re_eigs > threshold;

if isvector(rotorSpeed) && ~isvector(eigenvalues)
    if size(eigenvalues, 2) == length(rotorSpeed)
        speed_matrix = repmat(rotorSpeed(:)', size(eigenvalues, 1), 1);
    elseif size(eigenvalues, 1) == length(rotorSpeed)
        speed_matrix = repmat(rotorSpeed(:), 1, size(eigenvalues, 2));
    else
        error('Dimensions of rotorSpeed and eigenvalues are inconsistent.');
    end
else
    speed_matrix = rotorSpeed;
end


%plot
figure('Name', 'Damping vs Rotor Speed', 'NumberTitle', 'off')
clf
hold on 

if any(stable_idx(:))
    plot(speed_matrix(stable_idx), log_dec(stable_idx), 'bo', 'MarkerSize', 4)
end

if any(unstable_idx(:))
    plot(speed_matrix(unstable_idx), log_dec(unstable_idx), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r')
end

h_stable = plot(nan, nan, 'bo', 'MarkerSize', 6, 'DisplayName', 'Stable');
h_unstable = plot(nan, nan, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Unstable');

hold off

xlabel('Rotor Speed [rpm]')
ylabel('Logarithmic Decrement \delta') 
title('Internal Damping vs Rotor Speed')
grid on

yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); 


legend([h_unstable, h_stable], 'Location', 'best');

valid_data = log_dec(isfinite(log_dec));
if ~isempty(valid_data)
    y_span = max(abs(valid_data(abs(valid_data) < 10))); 
    if isempty(y_span) || y_span == 0, y_span = 1; end
    ylim([-0.5, max(1.5, y_span * 1.1)]); 
end

end

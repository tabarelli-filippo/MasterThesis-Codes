function [] = plotLogDecrement(rotorSpeed, eigenvalues, options)
%PLOTDAMPINGVSSPEED plots: - Logarithmic Decrement vs Rotor Speed
%
%INPUT: rotorSpeed      Vector of rotor speeds [RPM]
%       eigenvalues     Matrix or Vector of eigenvalues (complex)
%       threshold_stab  Stability threshold (default 1e-6)

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
        error('Dimensioni di rotorSpeed e eigenvalues non consistenti.');
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
function [] = waterfallPlot(X_freq, Y_speed_Hz, Z_plot, NX)
% WATERFALLPLOT  Generates a 3D waterfall (cascade) diagram showing the
%   spectral amplitude of the rotor response as a function of both
%   frequency and rotor speed.
%
%   The waterfall diagram is the time-frequency representation most
%   commonly used in run-up/coast-down analysis. Each row of Z_plot
%   corresponds to the frequency spectrum at one rotor speed; the
%   rows are stacked along the speed axis to form the 3D cascade.
%   Engine-order lines (NX) are optionally projected onto the base
%   plane (z = lower_limit) as diagonal reference lines.
%
% SYNTAX
%   waterfallPlot(X_freq, Y_speed_Hz, Z_plot)
%   waterfallPlot(X_freq, Y_speed_Hz, Z_plot, NX)
%
% INPUT ARGUMENTS
%   X_freq     - (vector, double) Frequency axis [Hz], corresponding to
%                the columns of Z_plot.
%   Y_speed_Hz - (vector, double) Rotor speed axis [Hz = rev/s],
%                corresponding to the rows of Z_plot.
%   Z_plot     - (nspeed x nfreq double) Amplitude matrix (typically
%                log-magnitude, e.g. 20*log10(|FFT|) or log10(amplitude)).
%                Rows = speed steps, columns = frequency bins.
%   NX         - (vector, double, default 1) Array of engine orders to
%                overlay as diagonal reference lines. Each NX(k) line
%                satisfies f = NX(k) * speed_Hz. Set NX = 0 to suppress.
%
% OUTPUT
%   A 3D figure titled 'Waterfall Diagram' is created with:
%     - Colour-interpolated mesh lines along the speed direction (rows)
%     - Colormap 'parula' clipped to [max(Z)-8, max(Z)] for contrast
%     - NX engine-order lines projected on the base plane in grey shades
%       (darker = higher order), labelled in the legend
%   The initial view angle is (azimuth=30°, elevation=45°).
%
% NOTES
%   - The amplitude lower limit is set 8 units below the global maximum
%     (or the data minimum if that is higher), clipping low-amplitude noise.
%   - For a 2D equivalent, use imagesc or pcolor on Z_plot directly.
%
% EXAMPLE
%   % Build spectrogram from run-up data
%   speeds_Hz = linspace(10, 80, 60);   % [Hz]
%   for ii = 1:length(speeds_Hz)
%       [f, Y, ~] = DFT(time, disp(:,ii), 2048);
%       Z(ii,:) = log10(abs(Y(1:N_half)));
%   end
%   waterfallPlot(f(1:N_half), speeds_Hz, Z, [1, 2, 3]);
%
% SEE ALSO
%   DFT, runUp, timeSimulation, plotCampbell

arguments (Input)
    X_freq (:,:) double {mustBeVector}
    Y_speed_Hz (:,:) double {mustBeVector}
    Z_plot (:,:) double
    NX (:,:) double {mustBeVector} = 1;
end

figure('Name','Waterfall Plot','NumberTitle', 'off');

% waterfall mesh: row-wise lines only
h = mesh(X_freq, Y_speed_Hz, Z_plot); 
set(h, 'FaceColor', 'none');
set(h, 'EdgeColor', 'interp');
set(h, 'MeshStyle', 'row');
set(h, 'LineWidth', 1);

xlabel('Frequency [Hz]');
ylabel('Rotor Speed [Hz]');
zlabel('Log Amplitude');
title('Waterfall Diagram');
grid on;
colormap(parula); 

max_Z = max(Z_plot, [], 'all');
lower_limit = max_Z - 8;
if lower_limit < min(Z_plot, [], 'all')
    lower_limit = min(Z_plot, [], 'all');
end
clim([lower_limit, max_Z]); 
zlim([lower_limit, max_Z]); 

%% Engine-order reference lines (projected onto the base plane)
hold on; 
legendHandles = [];

z_line_val = lower_limit;

if NX ~= 0
    speed_vec = Y_speed_Hz'; 
    z_vec = z_line_val * ones(size(speed_vec));

    freqLines = NX(:) * speed_vec;     % f = NX * speed_Hz
    
    hNX = gobjects(length(NX), 1);
    numLines = length(NX);
    
    % grey shading: lighter for lower orders, darker for higher orders
    maxGray = 0.7;
    minGray = 0.1;
    grayLevels = linspace(maxGray, minGray, numLines);
    
    for kk = 1:length(NX)
        lineStyle = '--';
        lineColor = [grayLevels(kk), grayLevels(kk), grayLevels(kk)];
        legendLabel = sprintf('%.1fX', NX(kk));  
        
        hNX(kk) = plot3(freqLines(kk, :), speed_vec, z_vec, ...
            'LineStyle', lineStyle, 'Color', lineColor, ...
            'LineWidth', 1.5, 'DisplayName', legendLabel);
    end

    legendHandles = [legendHandles; hNX];
    legend(legendHandles, 'Location', 'eastoutside');
end

hold off;
view(30, 45);

end

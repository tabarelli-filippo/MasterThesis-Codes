function [] = waterfallPlot(X_freq, Y_speed_Hz, Z_plot, NX)
%WATERFALLPLOT Genera un diagramma a cascata (Campbell/Waterfall)
% INPUT:
%   X_freq:     Vettore delle frequenze [Hz] (corrisponde alle colonne di Z)
%   Y_speed_Hz: Vettore delle velocità [Hz] (corrisponde alle righe di Z)
%   Z_plot:     Matrice delle ampiezze (es. log magnitude)
%   NX:         (Opzionale) Vettore degli ordini motore da tracciare

arguments (Input)
    X_freq (:,:) double {mustBeVector}
    Y_speed_Hz (:,:) double {mustBeVector} % CORRETTO: "double" (era "duoble")
    Z_plot (:,:) double
    NX (:,:) double {mustBeVector} = 1;
end

figure('Name','Waterfall Plot','NumberTitle', 'off');

% waterfall lines
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

%% Intersection lines
hold on; 
legendHandles = [];

z_line_val = lower_limit; 

if NX ~= 0
    speed_vec = Y_speed_Hz'; 
    z_vec = z_line_val * ones(size(speed_vec));

    freqLines = NX(:) * speed_vec; 
    
    hNX = gobjects(length(NX), 1);
    numLines = length(NX);
    
    % Scala di grigi
    maxGray = 0.7;
    minGray = 0.1;
    grayLevels = linspace(maxGray, minGray, numLines);
    
    for kk = 1:length(NX)
        lineStyle = '--';
        lineColor = [grayLevels(kk), grayLevels(kk), grayLevels(kk)];
        legendLabel = sprintf('%.1fX', NX(kk));  
        
        hNX(kk) = plot3(freqLines(kk, :), speed_vec, z_vec, 'LineStyle', lineStyle, 'Color', lineColor, ...
            'LineWidth', 1.5, 'DisplayName', legendLabel);
    end

    legendHandles = [legendHandles; hNX];
    legend(legendHandles, 'Location', 'eastoutside');
end

hold off;
view(30, 45);

end
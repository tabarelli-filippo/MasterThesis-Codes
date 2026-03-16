function [] = plotFRF(rotorSpeed,response,outnode)
%PLOTFRF plot the frequency response of the rotor system. Draws a figure
%for each (node,direction) specified in the input system
%
%INPUT: rotorSpeed  Array of speeds
%       response    Matrix with the response
%       outnode     Node and direction. EX: 1.2 = (Node 1, Dir y)
%                   Admissible directions are (1,2,3,4), respectively (x,y,theta,phi)

arguments (Input)
    rotorSpeed (:,:) double {mustBeVector}
    response (:,:) double
    outnode (:,:) double
end

if isscalar(rotorSpeed)
    error('error in plotFRF.m - there must be more than one rotor spin speed')
end

% admissible directions
adm_dirs = [1, 2, 3, 4];
nodes = fix(outnode);
dirs = round(rem(outnode, 1)*10);
if any(~ismember(dirs,adm_dirs),'all')
    error('error in plotFRF.m - directions must be admissible')
end

n_elements = length(outnode);
for ii = 1:n_elements
    
    idx = 4 * nodes(ii) - 4 + dirs(ii);
    responseData = response(idx, :);
    magnitudeData = abs(responseData);

    figure('Name', 'FRF (Bode Plot)', 'NumberTitle', 'off', 'Color', 'w');

    ax1 = subplot(2, 1, 1);
    plot(rotorSpeed, magnitudeData, 'r-', 'LineWidth', 1.5);
    set(gca, 'YScale', 'log');
    title('Magnitude Response');
    xlabel('Rotor Speed [rpm]');
    ylabel('Magnitude [m] (Log scale)');
    grid on;
    grid minor;

    y_max = max(magnitudeData);
    y_min_nonzero = min(magnitudeData(magnitudeData > 0));

    if ~isempty(y_min_nonzero) && ~isempty(y_max) && y_max > y_min_nonzero
        ylim([y_min_nonzero * 0.5, y_max * 1.5]);
    else
        ylim([1e-6, 1]);
    end

    ax2 = subplot(2, 1, 2);
    phaseData_deg = 180/pi * angle(responseData);
    plot(rotorSpeed, phaseData_deg, 'b-', 'LineWidth', 1.5);
    hold on;
    lineColor = [0.7 0.7 0.7];
    lineStyle = '--';
    yline([180, -180], 'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', 1);
    hold off;
    title('Phase response');
    xlabel('Rotor Speed [rpm]');
    ylabel('Phase [°]');
    grid on;
    grid minor;

    linkaxes([ax1, ax2], 'x');
    dirsAxis = {'\itx','\ity','\theta','\psi'};
    sgtitle(['Node Response: ', num2str(nodes(ii)), ', Direction: ', dirsAxis{dirs(ii)}]);
end
end
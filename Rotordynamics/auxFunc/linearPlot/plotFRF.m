function [] = plotFRF(rotorSpeed, response, outnode)
% PLOTFRF  Plots the Frequency Response Function (Bode plot) of the rotor
%   system as a function of rotational speed, for specified output nodes
%   and directions.
%
%   For each (node, direction) pair specified in outnode, a separate figure
%   is created showing the response amplitude (log scale) and phase versus
%   rotor speed. The amplitude is plotted in metres for translational DOFs
%   and radians for rotational DOFs.
%
% SYNTAX
%   plotFRF(rotorSpeed, response, outnode)
%
% INPUT ARGUMENTS
%   rotorSpeed - (vector, double) Array of rotor speeds [rpm] at which the
%                response was evaluated. Must have more than one element.
%   response   - (ndof x nspeed complex double) Complex displacement
%                response matrix as returned by FRF. Each column is the
%                response at one speed.
%   outnode    - (vector, double) Encoded list of output (node, direction)
%                pairs. Each entry is formatted as:
%                    node + direction/10
%                For example: 2.1 = node 2, x-direction
%                             3.2 = node 3, y-direction
%                Admissible directions:
%                    1 → x-displacement  (u)
%                    2 → y-displacement  (v)
%                    3 → x-rotation      (theta_u)
%                    4 → y-rotation      (theta_v)
%
% OUTPUT
%   One figure per (node, direction) pair, each containing:
%     Upper panel : amplitude |response| vs rotor speed [rpm], log y-scale
%     Lower panel : phase angle(response) in degrees vs rotor speed [rpm]
%   X-axes of both panels are linked. The figure title reports the node
%   number and direction symbol.
%
% NOTES
%   - The DOF index is computed as: idx = 4*(node-1) + direction
%   - Phase is wrapped to [-180°, 180°]; ±180° reference lines are shown.
%   - Amplitude axis limits are set automatically based on the data range,
%     with a fallback of [1e-6, 1] for near-zero responses.
%
% EXAMPLE
%   speeds_rpm = linspace(500, 5000, 500);
%   speeds_rads = speeds_rpm * pi/30;
%   resp = FRF(Rotor, speeds_rads);
%   % Plot x-response at node 3 and y-response at node 5
%   plotFRF(speeds_rpm, resp, [3.1, 5.2]);
%
% SEE ALSO
%   FRF, plotDisplacement, plotOrbit, critSpeeds

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

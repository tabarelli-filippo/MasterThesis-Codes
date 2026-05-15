function [] = plotOrbit(Mode, outnodes, eigenvalue, options)
% PLOTORBIT  Plots the precession orbits at specified rotor nodes for a
%   given mode shape and its associated eigenvalue.
%
%   Each orbit is the closed curve traced by the shaft centre at a node
%   during one precession cycle. It is obtained by evaluating the complex
%   displacement:
%       [u(alpha); v(alpha)] = Re(Mode(dof) * exp(i*alpha)),  alpha in [0, 2*pi)
%
%   The orbit is an ellipse whose shape (eccentricity and orientation)
%   encodes the whirl type: a circle indicates pure forward or backward
%   whirl; a line indicates a planar (degenerate) motion.
%
%   The traversal direction is indicated by:
%     x marker  : starting point (alpha = 0)
%     diamond   : ending point (alpha ≈ 2*pi)
%
%   Multiple nodes are superimposed on the same plot using different line
%   styles (solid, dashed, dotted, dash-dot) cycling through nodes.
%
% SYNTAX
%   plotOrbit(Mode, outnodes, eigenvalue)
%   plotOrbit(Mode, outnodes, eigenvalue, "titletext", txt)
%
% INPUT ARGUMENTS
%   Mode       - (ndof x 1 complex double) Mode shape eigenvector as
%                returned by charRoots (physical partition). DOF ordering:
%                [u1, v1, theta_u1, theta_v1, u2, ...]
%   outnodes   - (vector, positive integers) List of node indices to plot.
%                Each node contributes one orbit curve to the figure.
%   eigenvalue - (complex scalar) Eigenvalue associated with the mode.
%                The sign of Im(eigenvalue) determines the traversal
%                direction of the orbit:
%                  Im < 0 → standard complex notation (FW or BW per kappa)
%                  Im > 0 → conjugate representation; sign of j is flipped
%
% NAME-VALUE OPTIONS
%   "titletext" - (string, default '') Additional descriptive text for the
%                 figure subtitle (e.g., mode number, frequency, speed).
%
% OUTPUT
%   A figure with a square axis showing the overlaid orbits in the (u, v)
%   plane. Amplitudes are scaled to micrometres [µm] (factor 1e6 applied).
%   The legend identifies each orbit by node number. Grid and equal axes
%   are applied for correct ellipse visualisation.
%
% EXAMPLE
%   speeds = 1000 * pi/30;   % 1000 rpm in rad/s
%   [evals, evecs] = charRoots(Rotor, speeds);
%   figure;
%   plotOrbit(evecs(:,1), [2, 4, 6], evals(1), ...
%             "titletext", "Mode 1 at 1000 rpm");
%
% SEE ALSO
%   charRoots, modeWhirl, plotMode, plotDisplacement, poincareMap

arguments (Input)
    Mode (:,:) double {mustBeVector}
    outnodes (:,:) double {mustBeVector}
    eigenvalue (1,1)
    options.titletext (1,1) string = '';
end
titletext = options.titletext;
j = 1i;

N_out = length(outnodes);
outputdofx = 4*outnodes-3*ones(1,N_out);
outputdofy = 4*outnodes-2*ones(1,N_out);
outputdof = sort([outputdofx outputdofy]);

alpha = 0:0.01:1.9*pi;

% flip imaginary unit if eigenvalue has positive imaginary part
if imag(eigenvalue) < 0
        j = -j;
end

rr = 1e6*real(Mode(outputdof)*exp(j*alpha));
rmax = max(max(abs(rr)));

%% plot orbits
plotHandles = gobjects(N_out, 1); 
legendEntries = cell(N_out, 1);
lineStyles = {'-', '--', ':', '-.'};
numStyles = length(lineStyles);

hold on
for i=1:N_out
    styleIdx = mod(i - 1, numStyles) + 1;
    currentStyle = lineStyles{styleIdx};

    i1 = 2*i-1;
    i2 = i1+1;
    h = plot(rr(i1,:), rr(i2,:),'LineStyle', currentStyle, 'Color', 'k');
    plot(rr(i1,1),rr(i2,1),'kx')       % start marker
    plot(rr(i1,end),rr(i2,end),'kd')   % end marker
    plotHandles(i) = h;
    legendEntries{i} = ['Orbit node', num2str(outnodes(i))];
end
sgtitle('Node Orbits (u-v)')
title(titletext);
description = 'Cross: start - Diamond: end';
subtitle(description);

grid on
axis square
set(gca, 'Color', 'none');
axis([-rmax rmax -rmax rmax])
xlabel('u [\mum]');
ylabel('v [\mum]');
legend(plotHandles, legendEntries, 'Location', 'bestoutside');
hold off
end

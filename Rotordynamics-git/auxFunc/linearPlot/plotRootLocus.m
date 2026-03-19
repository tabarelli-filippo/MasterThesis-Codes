function [] = plotRootLocus(rotorSpeed, eigenvalues, options)
% PLOTROOTLOCUS  Plots the root locus and stability map of a rotor-bearing
%   system as functions of rotational speed.
%
%   Two subplots are generated:
%     1. Re(lambda) vs rotor speed: shows how the growth/decay rate of each
%        mode evolves with speed. Modes crossing Re = 0 become unstable.
%     2. Classic root locus: Im(lambda) vs Re(lambda), colour-coded by
%        speed. Each trajectory starts at low speed and ends at high speed;
%        trajectories crossing the imaginary axis indicate instability onset.
%
%   In both plots, stable modes (Re <= threshold) are shown with open
%   markers and unstable modes (Re > threshold) with filled markers.
%
% SYNTAX
%   plotRootLocus(rotorSpeed, eigenvalues)
%   plotRootLocus(rotorSpeed, eigenvalues, "threshold_stab", thr)
%
% INPUT ARGUMENTS
%   rotorSpeed  - (vector, double) Rotor speeds at which eigenvalues were
%                 evaluated [rpm]. Size must be consistent with eigenvalues.
%   eigenvalues - (matrix, complex double) Complex eigenvalues as returned
%                 by charRoots. Accepts both (ncdof x nspeed) and flattened
%                 layouts; the function uses element-wise Re/Im extraction.
%
% NAME-VALUE OPTIONS
%   "threshold_stab" - (double, default 1e-6) Stability threshold on
%                      Re(lambda) [1/s]. Eigenvalues with Re > threshold
%                      are classified as unstable.
%
% OUTPUT
%   A figure titled 'Root Locus and Stability' with two subplots:
%     Upper : Re(lambda) [1/s] vs rotor speed [rpm]
%             Blue circles = stable; red filled circles = unstable;
%             dashed black line = Re = 0 boundary.
%     Lower : Im(lambda) [rad/s] vs Re(lambda) [1/s], coloured by speed
%             using the 'parula' colormap; a colorbar shows the speed scale.
%             Dashed black line = imaginary axis (Re = 0).
%
% EXAMPLE
%   speeds_rpm = linspace(500, 8000, 200);
%   speeds_rads = speeds_rpm * pi/30;
%   [evals, ~] = charRoots(Rotor, speeds_rads);
%   plotRootLocus(speeds_rpm, evals);
%
% SEE ALSO
%   charRoots, plotLogDecrement, plotCampbell, LogDecrement

arguments (Input)
    rotorSpeed (:,:) double {mustBeVector}
    eigenvalues (:,:) double
    options.threshold_stab (1,1) double = 1e-6
end
threshold = options.threshold_stab;
figure('Name', 'Root Locus and Stability','NumberTitle', 'off')
clf

re_eigs = real(eigenvalues);
im_eigs = imag(eigenvalues);

stable_idx = re_eigs <= threshold;
unstable_idx = re_eigs > threshold;

y_max = max(max(abs(re_eigs)), 1e-3);

%% Subplot 1: Re(lambda) vs rotor speed — stability map
subplot(2, 1, 1)
hold on 

plot(rotorSpeed(stable_idx), re_eigs(stable_idx), 'bo', 'MarkerSize', 4)
plot(rotorSpeed(unstable_idx), re_eigs(unstable_idx), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r')

hold off
xlabel('Rotor Speed [rpm]')
ylabel('Re(\lambda) [1/s]')
title('Stability vs Rotor Speed')
grid on

ax = gca;
line(ax.XLim, [0 0], 'Color', 'k', 'LineStyle', '--')
ylim([-y_max, y_max]);

%% Subplot 2: root locus — Im vs Re, coloured by speed
subplot(2, 1, 2)
hold on 

h1 = scatter(re_eigs(stable_idx),   im_eigs(stable_idx),   25, rotorSpeed(stable_idx)); 
h2 = scatter(re_eigs(unstable_idx), im_eigs(unstable_idx), 25, rotorSpeed(unstable_idx), 'filled'); 
colormap('parula');
c = colorbar;
c.Label.String = 'Rotor Speed [rpm]';

hold off
xlabel('Re(\lambda) [1/s]')
ylabel('Im(\lambda) [rad/s]')
title('Root Locus')
grid on

ax = gca;
line([0 0], ax.YLim, 'Color', 'k', 'LineStyle', '--')

xlim([-y_max, y_max]);
handle = [h1; h2];
leg_text = {'Empty = Stable','Filled = Unstable'};
legend(handle,leg_text);
end

function [] = poincareMap(Time, x, y, rotorSpeed)
% POINCAREMAP  Plots the Poincaré map of a rotor orbit at a given
%   rotational speed, showing both the sampled section points and the
%   full continuous orbit for comparison.
%
%   The Poincaré section is defined by sampling the (x, y) orbit once per
%   shaft revolution, at time instants:
%       t_k = T_start + k * T,    T = 2*pi / Omega
%
%   The sampled points are obtained by linear interpolation of the
%   time-domain trajectory. The Poincaré map is a classical tool for
%   distinguishing between:
%     - Period-1 response   : a single fixed point in the map
%     - Period-N response   : N discrete points
%     - Quasi-periodic      : a closed curve (torus section)
%     - Chaotic response    : a fractal/scattered cloud
%
% SYNTAX
%   poincareMap(Time, x, y, rotorSpeed)
%
% INPUT ARGUMENTS
%   Time       - (ntime x 1 double) Time vector [s], monotonically
%                increasing. Typically from timeSimulation or runUp.
%   x          - (ntime x 1 double) x-displacement time series [m]
%                (e.g., u-direction DOF extracted from the q output)
%   y          - (ntime x 1 double) y-displacement time series [m]
%                (e.g., v-direction DOF extracted from the q output)
%   rotorSpeed - (scalar double) Constant rotor speed [rad/s] used to
%                compute the sampling period T = 2*pi / rotorSpeed
%
% OUTPUT
%   A figure titled 'Poincaré Map Analysis' with two side-by-side subplots:
%     Left panel  : Poincaré section — sampled (x, y) points, coloured
%                   by sampling time using the default colormap. A cross
%                   at the origin marks the shaft reference position.
%                   A colorbar shows the time axis [s].
%     Right panel : Full orbit trajectory (thin black line) overlaid with
%                   the Poincaré sample points (red dots). A cross at
%                   the origin is also shown.
%   Both panels use equal aspect ratio and display rotor speed in the title.
%
% EXAMPLE
%   [t, q, ~] = timeSimulation(Rotor, [0, 5], 500, zeros(ndof,1), zeros(ndof,1));
%   node = 3;
%   x_sig = q(:, 4*node-3);   % x-displacement at node 3
%   y_sig = q(:, 4*node-2);   % y-displacement at node 3
%   poincareMap(t, x_sig, y_sig, 500);
%
% SEE ALSO
%   timeSimulation, runUp, plotBifurcation, plotOrbit, DFT

arguments (Input)
    Time (:,1) double {mustBeVector}
    x (:,1) double {mustBeVector}
    y (:,1) double {mustBeVector}
    rotorSpeed (1,1) double
end

T = 2*pi / rotorSpeed;                            % one revolution period [s]
t_poincare = (Time(1):T:Time(end))';              % Poincaré section times

x_p = interp1(Time, x, t_poincare, 'linear', 'extrap');
y_p = interp1(Time, y, t_poincare, 'linear', 'extrap');

figure('Name', 'Poincaré Map Analysis', 'NumberTitle', 'off');

%% Left panel: Poincaré section coloured by time
subplot(1,2,1);
scatter(x_p, y_p, 10, t_poincare, 'filled');
hold on;
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
hold off 
xlabel('u [m]');
ylabel('v [m]');
title(sprintf('Poincaré Map (Speed: %.0f rpm)', rotorSpeed*60/2/pi));
pbaspect([1 1 1]);
grid on; grid minor; box on;
colorbar; ylabel(colorbar, 'Time [s]');

%% Right panel: full orbit with Poincaré samples overlaid
subplot(1,2,2);
hold on;
plot(x, y, 'k', 'LineWidth', 0.5);                     % full orbit
plot(x_p, y_p, 'r.', 'MarkerSize', 10);                % Poincaré points
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);  % shaft axis
hold off
xlabel('u [m]');
ylabel('v [m]');
title('Full Orbit and Samples Location');
pbaspect([1 1 1]);
grid on; grid minor; box on;

end

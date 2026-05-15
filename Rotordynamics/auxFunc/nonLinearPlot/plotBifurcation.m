function [] = plotBifurcation(rotorSpeed, Time, dispData)
% PLOTBIFURCATION  Generates a bifurcation diagram for a rotor system by
%   sampling the displacement response at each shaft revolution (Poincaré
%   section) and plotting the sampled values against rotor speed.
%
%   For each speed in rotorSpeed, the Poincaré section is defined as the
%   set of times spaced exactly one shaft revolution period apart:
%       T = 2*pi / Omega,   t_k = T_start + k*T
%
%   The displacement at each Poincaré instant is obtained by linear
%   interpolation of the time-domain signal. The resulting (speed, disp)
%   pairs are then plotted as a dot cloud.
%
%   Bifurcation diagrams reveal qualitative changes in the attractor
%   structure as speed varies:
%     - A single point per speed   → period-1 (synchronous) response
%     - N discrete points          → period-N subharmonic response
%     - A dense cloud              → quasi-periodic or chaotic response
%
% SYNTAX
%   plotBifurcation(rotorSpeed, Time, dispData)
%
% INPUT ARGUMENTS
%   rotorSpeed - (1 x nspeed double) Array of rotor speeds [rad/s] at
%                which the time-domain response was computed. Each entry
%                corresponds to a column of dispData.
%   Time       - (ntime x 1 double) Uniform (or non-uniform) time vector
%                [s] shared by all speed simulations. Must be monotonically
%                increasing.
%   dispData   - (ntime x nspeed double) Displacement time series at each
%                speed. Typically a single DOF extracted from the output of
%                timeSimulation or runUp.
%
% OUTPUT
%   A figure titled 'Bifurcation Analysis' is created showing a dot plot
%   of Poincaré-sampled displacement [m] versus rotor speed [rad/s].
%   X-axis limits are set to ±5% beyond the speed range.
%
% NOTES
%   - Speeds for which no valid Poincaré points fall within [Time(1), Time(end)]
%     are silently skipped.
%   - NaN interpolation results (caused by extrapolation gaps) are removed.
%   - For a single speed, use poincareMap instead for a more detailed view.
%
% EXAMPLE
%   speeds = linspace(100, 600, 80);   % [rad/s]
%   [t, q, ~] = timeSimulation(Rotor, [0, 2], speeds(1), zeros(ndof,1), zeros(ndof,1));
%   % (collect dispData across speeds in a loop, then:)
%   plotBifurcation(speeds, t, dispData);
%
% SEE ALSO
%   poincareMap, timeSimulation, runUp, DFT

arguments (Input)
    rotorSpeed (1,:) double {mustBeVector} 
    Time (:,1) double {mustBeVector}       
    dispData (:,:) double                  
end

nspeed = length(rotorSpeed);
X_plot = [];
Y_plot = [];

for ii = 1:nspeed
    omega = rotorSpeed(ii);
    T = 2*pi / omega;           % one shaft revolution period
    
    t_poincare = Time(1):T:Time(end);
    
    if isempty(t_poincare)
        continue; 
    end
    
    y_current = dispData(:, ii);
    points = interp1(Time, y_current, t_poincare, 'linear');
    
    valid_mask = ~isnan(points);
    points = points(valid_mask);
    
    x_vals = repmat(omega, length(points), 1);
    
    X_plot = [X_plot; x_vals];
    Y_plot = [Y_plot; points(:)]; 
end

figure('Name','Bifurcation Analysis','NumberTitle', 'off', 'Color', 'w');
plot(X_plot, Y_plot, 'k.', 'MarkerSize', 3); 
grid on;
grid minor;
xlabel('Rotor Speed [rad/s]', 'FontWeight', 'bold');
ylabel('Poincaré Displacement [m]', 'FontWeight', 'bold');
title('Bifurcation Plot');
xlim([min(X_plot)*0.95, max(X_plot)*1.05]);

end

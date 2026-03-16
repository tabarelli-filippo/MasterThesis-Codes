function [] = poincareMap(Time, x, y, rotorSpeed)
%POINCAREMAP: plots Poincaré map for a given rotorSpeed
arguments (Input)
    Time (:,1) double {mustBeVector}
    x (:,1) double {mustBeVector}
    y (:,1) double {mustBeVector}
    rotorSpeed (1,1) double % [rad/s]
end

T = 2*pi / rotorSpeed;
t_poincare = (Time(1):T:Time(end))';

x_p = interp1(Time, x, t_poincare, 'linear', 'extrap');
y_p = interp1(Time, y, t_poincare, 'linear', 'extrap');

figure('Name', 'Poincaré Map Analysis', 'NumberTitle', 'off');

%% subplot 1
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

%% subplot 2
subplot(1,2,2);
hold on;
plot(x, y, 'k', 'LineWidth', 0.5);
plot(x_p, y_p, 'r.', 'MarkerSize', 10);
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
hold off

xlabel('u [m]');
ylabel('v [m]');

title('Full Orbit and Samples Location');
pbaspect([1 1 1]);
grid on; grid minor; box on;

end
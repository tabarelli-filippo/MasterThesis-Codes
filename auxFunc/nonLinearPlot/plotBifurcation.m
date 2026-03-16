function [] = plotBifurcation(rotorSpeed, Time, dispData)
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
    T = 2*pi / omega;
    
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
xlabel('Rotor Speed [RPM]', 'FontWeight', 'bold');
ylabel('Poincaré Displacement [m]', 'FontWeight', 'bold');
title('Bifurcation Plot');
xlim([min(X_plot)*0.95, max(X_plot)*1.05]);

end
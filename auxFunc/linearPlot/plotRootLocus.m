function [] = plotRootLocus(rotorSpeed,eigenvalues,options)
%PLOTROOTLOCUS plots:   -Real(eigs) vs Speed to evaluate stability over rotor
%                       speed
%                       -Root Locus: Re(eigs) vs Imag(eigs)
%
%INPUT: rotorSpeed      Vector of rotor speeds
%       eigenvalues     List of eigenvalues 
%       threshold_stab  stability threshold caused by numerical errors

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

%% Real(eigs) vs Speed
subplot(2, 1, 1)
hold on 

plot(rotorSpeed(stable_idx), re_eigs(stable_idx), 'bo', 'MarkerSize', 4)
plot(rotorSpeed(unstable_idx), re_eigs(unstable_idx), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r')

hold off
xlabel('Rotor Speed [rpm]')
ylabel('Re(eigs) ')
title('Stability vs Rotor Speed')
grid on

ax = gca;
line(ax.XLim, [0 0], 'Color', 'k', 'LineStyle', '--')
ylim([-y_max, y_max]);

%% Root Locus
subplot(2, 1, 2)
hold on 

h1 = scatter(re_eigs(stable_idx), im_eigs(stable_idx), 25, rotorSpeed(stable_idx)); 
h2 = scatter(re_eigs(unstable_idx), im_eigs(unstable_idx), 25, rotorSpeed(unstable_idx), 'filled'); 
colormap('parula');
c = colorbar;
c.Label.String = 'Rotor Speed [RPM]';

hold off
xlabel('Real(eigs)')
ylabel('Imag(eigs)')
title('Root Locus')
grid on

ax = gca;
line([0 0], ax.YLim, 'Color', 'k', 'LineStyle', '--')

xlim([-y_max, y_max]);
handle = [h1; h2];
leg_text = {'Empty = Stable','Filled = Unstable'};
legend(handle,leg_text);
end
function [] = plotOrbit(Mode,outnodes,eigenvalue,options)
%PLOTORBIT plots di orbit of a given node list. Cross is the starting
%point, diamond is the ending point.
%
%INPUT: Mode        Vector with the modes associated to the eigenvalue
%       outnodes    List of nodes to plot
%       eigenvalue  eigenvalue associated to the mode
%       titletext   auxiliary text description
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

if imag(eigenvalue)<0
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
    plot(rr(i1,1),rr(i2,1),'kx')
    plot(rr(i1,end),rr(i2,end),'kd')
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
legend(plotHandles, legendEntries, 'Location', 'bestoutside');
hold off
end
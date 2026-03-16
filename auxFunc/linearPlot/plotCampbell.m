function [] = plotCampbell(RotorSpeed,eigvals,NX,isDamped,kappa)
%PLOTCAMPBELL plots the Campbell Diagram given range of speeds
%INPUT: RotorSpeed  vector of rotor spin speeds [rad/s] - must be more
%                   than one speed
%
%       eigvals     complex eigenvalues, obtained from charRoot.m
%
%       NX          array of armonic intersections. NX=0 plots no line
%
%       isDamped    parameter to plot undamped natural frequencies (default,
%                   dampedFreq=0) or undamped natural frequencies (dampedFreq=1)
%
%       kappa       rotor orbit and direction, obtained from charRoot.m
%                   k>0 FW mode, k>0 BW mode

arguments (Input)
    RotorSpeed (:,:) double {mustBeVector}
    eigvals (:,:) double
    NX (1,:) double = 0 %default
    isDamped (1,1) logical = false %default
    kappa (:,:,:) double = []
end

RotorSpeed_rpm = RotorSpeed*30/pi;
figure('Name','Campbell Diagram','NumberTitle', 'off')
clf
hold on; 
grid on; grid minor;
xlabel('Rotor Speed [rpm]');

% y-axis damped natural frequencies or natural frequencies
if isDamped
    natFreqs_Hz = abs(imag(eigvals))/(2*pi);  % damped natural frequencies
    ylabel('Damped Frequency [Hz]');
else
    natFreqs_Hz = abs(eigvals)/(2*pi);  % natural frequencies
    ylabel('Natural Frequency [Hz]');
end
max_Freq = max(natFreqs_Hz,[],'all');
ylim([0, 1.1*max_Freq])

title('Campbell Diagram');
if isempty(kappa)
    plot(RotorSpeed_rpm, natFreqs_Hz, 'k')
    legendHandles= [];
else
    [~,n2,nspeed] = size(kappa);
    dotFW  = nan(n2, nspeed);
    dotBW  = nan(n2, nspeed);
    dotMix = nan(n2, nspeed);

    % All k => 0 : forward mode
    maskFW = squeeze(all(kappa >= 0, 1));

    % All k < 0 : backward mode
    maskBW = squeeze(all(kappa <= 0, 1)) & ~maskFW;

    % mixed mode
    maskMix = ~maskFW & ~maskBW;

    dotFW(maskFW)   = natFreqs_Hz(maskFW);
    dotBW(maskBW)   = natFreqs_Hz(maskBW);
    dotMix(maskMix) = natFreqs_Hz(maskMix);
    plot(RotorSpeed_rpm, natFreqs_Hz, 'k');
    plot(RotorSpeed_rpm, dotFW, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(RotorSpeed_rpm, dotBW, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(RotorSpeed_rpm, dotMix, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
    
    % dummies for legend
    hFW  = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'DisplayName', 'Forward Modes');
    hBW  = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Backward Modes');
    hMix = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Mixed Modes');
    
    legendHandles = [hFW, hBW, hMix];
end

% intersection lines

if NX ~= 0
    freqLines = NX'*RotorSpeed/(2*pi);
    hNX = gobjects(length(NX), 1);
    
    numLines = length(NX);
    maxGray = 0.7;
    minGray = 0.1;
    grayLevels = linspace(maxGray, minGray, numLines);

    for kk = 1:length(NX)
        lineStyle = '--';
        lineColor = [grayLevels(kk), grayLevels(kk), grayLevels(kk)];
        legendLabel = sprintf('%.2fX', NX(kk));  
        hNX(kk) = plot(RotorSpeed_rpm, freqLines(kk, :),'LineStyle', lineStyle,'Color', lineColor, 'DisplayName', legendLabel);
    end
    legendHandles = [legendHandles, hNX'];
end

legend(legendHandles, 'Location', 'eastoutside');
hold off

end


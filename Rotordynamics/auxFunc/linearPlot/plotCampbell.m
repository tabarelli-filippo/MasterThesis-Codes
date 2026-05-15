function [] = plotCampbell(RotorSpeed,eigvals,NX,isDamped,kappa)
% PLOTCAMPBELL  Generates the Campbell diagram (interference diagram) of a
%   rotor-bearing system over a range of rotational speeds.
%
%   The Campbell diagram plots natural (or damped) frequencies as functions
%   of rotor speed, overlaid with engine-order excitation lines (NX lines).
%   Intersections between natural frequency curves and NX lines identify
%   potential critical speeds. When the whirl ratio kappa is supplied, each
%   data point is colour-coded to indicate the direction of precession.
%
% SYNTAX
%   plotCampbell(RotorSpeed, eigvals)
%   plotCampbell(RotorSpeed, eigvals, NX)
%   plotCampbell(RotorSpeed, eigvals, NX, isDamped)
%   plotCampbell(RotorSpeed, eigvals, NX, isDamped, kappa)
%
% INPUT ARGUMENTS
%   RotorSpeed - (vector, double) Rotor speeds [rad/s] at which eigenvalues
%                were computed. Must have more than one element.
%   eigvals    - (ncdof x nspeed complex double) Complex eigenvalues as
%                returned by charRoots. Natural frequencies are extracted
%                as |Im(lambda)| (damped) or |lambda| (undamped).
%   NX         - (1 x nX double, default 0) Array of engine orders for
%                excitation lines. Each value NX(k) plots the line
%                f = NX(k) * Omega / (2*pi). Set NX = 0 to suppress lines.
%   isDamped   - (logical, default false) false → plot undamped natural
%                frequencies |lambda|; true → plot damped natural
%                frequencies |Im(lambda)|.
%   kappa      - (ndof x ncdof x nspeed double, optional) Whirl ratio
%                tensor from charRoots. When supplied, mode markers are
%                colour-coded:
%                  Blue  = Forward Whirl (all kappa >= 0 at the mode's nodes)
%                  Red   = Backward Whirl (all kappa <= 0)
%                  Black = Mixed (some nodes FW, some BW)
%
% OUTPUT
%   A figure titled 'Campbell Diagram' is created. NX excitation lines are
%   drawn in shades of grey (darkest = highest order) and labelled in the
%   legend. The legend is placed outside the axes on the right.
%
% EXAMPLE
%   speeds = linspace(100, 2000, 150);      % [rad/s]
%   [evals, ~, kap] = charRoots(Rotor, speeds);
%   [evals, ~, kap] = sortModesMAC(evals, [], kap);
%   plotCampbell(speeds, evals, [1, 2, 3], false, kap);
%
% SEE ALSO
%   charRoots, sortModesMAC, critSpeeds, plotLogDecrement, plotRootLocus

arguments (Input)
    RotorSpeed (:,:) double {mustBeVector}
    eigvals (:,:) double
    NX (1,:) double = 0
    isDamped (1,1) logical = false
    kappa (:,:,:) double = []
end

RotorSpeed_rpm = RotorSpeed*30/pi;
figure('Name','Campbell Diagram','NumberTitle', 'off')
clf
hold on; 
grid on; grid minor;
xlabel('Rotor Speed [rpm]');

% y-axis: damped or undamped natural frequencies
if isDamped
    natFreqs_Hz = abs(imag(eigvals))/(2*pi);
    ylabel('Damped Frequency [Hz]');
else
    natFreqs_Hz = abs(eigvals)/(2*pi);
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

    % Forward whirl: all nodal kappa >= 0
    maskFW = squeeze(all(kappa >= 0, 1));
    % Backward whirl: all nodal kappa <= 0 (and not all zero / FW)
    maskBW = squeeze(all(kappa <= 0, 1)) & ~maskFW;
    % Mixed: neither purely FW nor purely BW
    maskMix = ~maskFW & ~maskBW;

    dotFW(maskFW)   = natFreqs_Hz(maskFW);
    dotBW(maskBW)   = natFreqs_Hz(maskBW);
    dotMix(maskMix) = natFreqs_Hz(maskMix);
    plot(RotorSpeed_rpm, natFreqs_Hz, 'k');
    plot(RotorSpeed_rpm, dotFW, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(RotorSpeed_rpm, dotBW, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(RotorSpeed_rpm, dotMix, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
    
    % dummy handles for legend
    hFW  = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'DisplayName', 'Forward Modes');
    hBW  = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Backward Modes');
    hMix = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Mixed Modes');
    
    legendHandles = [hFW, hBW, hMix];
end

% engine-order excitation lines
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

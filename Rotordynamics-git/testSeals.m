% TESTSEALS  Test and visualisation script for the labyrinth seal model.
%
%   Evaluates the rotordynamic coefficients (K, k, C, c) of a labyrinth
%   seal as a function of inlet swirl velocity for both choked and unchoked
%   flow conditions, using the Scharrer-Childs bulk-flow model implemented
%   in labySeals.
%
%   The seal geometry and operating conditions are taken from published
%   validation data (teeth-on-stator configuration). The script iterates
%   over a range of inlet swirl velocities U_inlet and plots all four
%   rotordynamic coefficients for both flow regimes on a single figure,
%   allowing direct comparison.
%
%   SEAL PARAMETERS (TOS, representative values):
%     Nt  = 16 teeth
%     Cr  = 0.406 mm  (radial clearance)
%     L   = 3.175 mm  (tooth pitch)
%     B   = 3.175 mm  (tooth depth)
%     Rs  = 76.5 mm   (seal radius)
%     T   = 300 K     (gas temperature)
%     nu  = 0.144e-4 m²/s (kinematic viscosity)
%     omega = 3000 rpm
%
%   FLOW CONDITIONS:
%     Choked:   P_in = 7.6 bar,  P_out = 0.943 bar
%     Unchoked: P_in = 2.9 bar,  P_out = 0.943 bar
%
%   FIGURE LAYOUT:
%     (1,1) Direct stiffness K [N/m] vs inlet swirl [m/s]
%     (1,2) Cross-coupled stiffness k [N/m] vs inlet swirl [m/s]
%     (2,1) Direct damping C [N·s/m] vs inlet swirl [m/s]
%     (2,2) Cross-coupled damping c [N·s/m] vs inlet swirl [m/s]
%
%   Blue = choked flow (P_in = 7.6 bar)
%   Red  = unchoked flow (P_in = 2.9 bar)
%
% SEE ALSO
%   labySeals, sealCoeffs, backwardPressureSolver, forwardPressureSolver

clear
close all
addpath(genpath('auxFunc/'))

%% Seal geometry and operating conditions
sealType = 'TOS';
P_out    = 0.943e5;     % downstream (sump) pressure [Pa]
Nt       = 16;          % number of teeth
H        = 0.000406;    % tooth radial clearance Cr [m]
L        = 0.003175;    % tooth pitch [m]
T        = 300;         % gas temperature [K]
Rs       = 0.0765;      % seal radius [m]
B        = 0.003175;    % tooth depth [m]
omega    = 3000*pi/30;  % shaft angular velocity [rad/s]
nu       = 0.144e-4;    % gas kinematic viscosity [m²/s]

%% Inlet swirl sweep
U_inlet = -40:10:40;    % tangential inlet velocity [m/s]
nn = length(U_inlet);
K = zeros(nn,1); k = zeros(nn,1);
C = zeros(nn,1); c = zeros(nn,1);

figure('Name', 'Rotordynamic Coefficients vs Inlet Swirl', 'Color', 'w', 'NumberTitle', 'off');
myColors = [0 0.447 0.741;    % blue  — choked
            0.85 0.325 0.098]; % orange-red — unchoked

for ii = 1:2
    switch ii
        case 1
            P_in = 7.6e5;       % choked inlet pressure [Pa]
            flowCondition = 'Choked';
        case 2
            P_in = 2.9e5;       % unchoked inlet pressure [Pa]
            flowCondition = 'Unchoked';
    end

    % compute rotordynamic coefficients at each inlet swirl value
    for jj = 1:nn
        [~, Kseal, Cseal] = labySeals(sealType, P_in, P_out, Nt, H, L, T, ...
            Rs, B, U_inlet(jj), omega, nu);
        K(jj) = Kseal(1,1);
        k(jj) = Kseal(1,2);
        C(jj) = Cseal(1,1);
        c(jj) = Cseal(1,2);
    end

    legendaStr = sprintf('P_{in} = %.2f bar [%s]', P_in/1e5, flowCondition);

    % --- Direct Stiffness K ---
    subplot(2,2,1); hold on;
    plot(U_inlet, K, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Direct Stiffness (K)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Stiffness [N/m]');
    xlim([min(U_inlet) max(U_inlet)]);

    % --- Cross-Coupled Stiffness k ---
    subplot(2,2,2); hold on;
    plot(U_inlet, k, '-s', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Cross-Coupled Stiffness (k)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Stiffness [N/m]');
    xlim([min(U_inlet) max(U_inlet)]);

    % --- Direct Damping C ---
    subplot(2,2,3); hold on;
    plot(U_inlet, C, '-^', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Direct Damping (C)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Damping [N·s/m]');
    xlim([min(U_inlet) max(U_inlet)]);

    % --- Cross-Coupled Damping c ---
    subplot(2,2,4); hold on;
    plot(U_inlet, c, '-d', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Cross-Coupled Damping (c)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Damping [N·s/m]');
    xlim([min(U_inlet) max(U_inlet)]);
end

for i = 1:4
    subplot(2,2,i);
    legend('show', 'Location', 'best');
    grid on; grid minor;
end

sgtitle('Influence of Inlet Swirl on Rotordynamic Coefficients');

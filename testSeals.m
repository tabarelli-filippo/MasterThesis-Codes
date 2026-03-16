clear
close all
addpath(genpath('auxFunc/'))

%% testing labySeals
sealType = 'TOS';
P_out = 0.943e5; %[Pa]
Nt = 16; % Number of teeth
H = 0.000406; % tooth heigth [m]
L = 0.003175; % Length [m]
T = 300; % Temperature [K]
Rs = 0.0765; % Seal radius [m]
B = 0.003175; % Width in meters
omega = 3000*pi/30; % [rad/s]
nu = 0.144e-4; % Kinematic viscosity [m^2/s]
U_inlet = -40:10:40;
nn = length(U_inlet);
K = zeros(nn,1);
k = zeros(nn,1);
C = zeros(nn,1);
c = zeros(nn,1);

figure('Name', 'Rotordynamic Coefficients vs Inlet Swirl', 'Color', 'w','NumberTitle', 'off');

myColors = [0 0.447 0.741; 0.85 0.325 0.098]; 

for ii=1:2
    switch ii
        case 1 % choked
            P_in = 7.6e5; %[Pa]
            flowCondition = 'Choked';
        case 2 % unchoked
            P_in = 2.9e5; %[Pa]
            flowCondition = 'Unchoked';
    end
    
    % Calcolo coefficienti
    for jj = 1:nn
        [mdot_leakage,Kseal,Cseal] = labySeals(sealType,P_in,P_out,Nt,H,L,T,Rs,B,U_inlet(jj),...
            omega,nu);
        K(jj) = Kseal(1,1);
        k(jj) = Kseal(1,2);
        C(jj) = Cseal(1,1);
        c(jj) = Cseal(1,2);
    end
    
    legendaStr = sprintf('P_{in} = %.2f bar [%s]', P_in/1e5, flowCondition);
    
    % --- PLOTTING ---
    
    % 1. Direct Stiffness (K)
    subplot(2,2,1)
    hold on; % Importante: mantiene il grafico precedente
    plot(U_inlet, K, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Direct Stiffness (K)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Stiffness [N/m]');
    xlim([min(U_inlet) max(U_inlet)]);

    % 2. Cross-Coupled Stiffness (k)
    subplot(2,2,2)
    hold on;
    plot(U_inlet, k, '-s', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Cross-Coupled Stiffness (k)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Stiffness [N/m]');
    xlim([min(U_inlet) max(U_inlet)]);

    % 3. Direct Damping (C)
    subplot(2,2,3)
    hold on;
    plot(U_inlet, C, '-^', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Direct Damping (C)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Damping [N s/m]');
    xlim([min(U_inlet) max(U_inlet)]);

    % 4. Cross-Coupled Damping (c)
    subplot(2,2,4)
    hold on;
    plot(U_inlet, c, '-d', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'Color', myColors(ii,:), 'DisplayName', legendaStr);
    title('Cross-Coupled Damping (c)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Inlet Tangential Velocity [m/s]');
    ylabel('Damping [N s/m]');
    xlim([min(U_inlet) max(U_inlet)]);
end

for i = 1:4
    subplot(2,2,i);
    legend('show', 'Location', 'best');
    grid on; grid minor;
end

sgtitle('Influence of Inlet Swirl on Rotordynamic Coefficients');
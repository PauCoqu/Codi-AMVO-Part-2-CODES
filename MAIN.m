clc;
clear all;
close all;

%% DEFINING THE PROBLEM

% Definition of variables for Wing
Nw = 512; % Number of spanwise segments for the wing
b = 12; %[m] Wingspan
cr = 1.8; %[m] Root chord length
ct = 1.2; %[m] Tip chord length
S_wing = (2*cr*(1+ct/cr+(ct/cr)^2))/(3*(1+ct/cr))*b; %[m^2] Wing surface

% Definition of variables for Canard
Nc = 512; % Number of spanwise segments for the canard
bh = 4.8; %[m] Canard span
crh = 0.9; %[m] Canard root chord length
cth = 0.6; %[m] Canard tip chord length
lh = 6; %[m] Longitudinal distance from wing to canard
S_canard = (2*crh*(1+cth/crh+(cth/crh)^2))/(3*(1+cth/crh))*bh; %[m^2] Canard surface

Q_inf = 1; % Velocity magnitude [m/s]
Cd_v = 0.0075;
Sv = 2.4; %[m^2]
iw = 0; %[rad] incidence for NACA2412 (wing)
ic = deg2rad(2); %[rad] incidence for NACA0010 (canard)
rho = 1; %[kg/m^3]

% Twist angles for study
theta_deg = -8:1:0;
theta = deg2rad(theta_deg);

% Load aerodynamic coefficients
[Cl_alpha_0010, Cl0_0010, Cm14_0010, Cl_alpha_2412, Cl0_2412, Cm14_2412, ...
 Cl_alpha_Fowler15, Cl0_Fowler15, Cm14_Fowler15] = AeroCoefficientsPart1();

% Generate wing and canard geometry
[X_wing, Xc_wing, c_wing, X_canard, Xc_canard, c_canard] = ...
    generateWingAndCanard(Nw, b, cr, ct, Nc, bh, crh, cth, lh, theta);

%% PART 1A - DEFINE ADEQUATE WING TWIST ANGLE
alpha = deg2rad(4); % Angle of attack
Ur = [-cos(alpha); 0; sin(alpha)];

% Initialize storage
Cl_w_iso_t = zeros(size(theta));
E_t_w = zeros(size(theta));
Cl_w_iso_dis = zeros(Nw, length(theta));
alpha_ind_w = zeros(Nw, length(theta));
Cd_visc_w = zeros(Nw, length(theta));
Cd_visc_t_w = zeros(1, length(theta));
Cd_induced_w = zeros(Nw, length(theta));

for i = 1:length(theta)
    t_wing = ComputeTwistDistribution(theta(i), Nw);
    gamma_w_isol = Gamma_wing(Nw, Q_inf, alpha, Ur, X_wing, Xc_wing, c_wing, iw, t_wing, Cl_alpha_2412, Cl0_2412);
    [Cl_w_iso_t(i), Cl_w_iso_dis(:,i), alpha_ind_w(:,i), ...
     Cd_visc_w(:,i), Cd_visc_t_w(i), Cd_induced_w(:,i), ...
     ~, E_t_w(i), ~, ~] = ...
     Lift_and_CL(Nw, gamma_w_isol, X_wing, Q_inf, S_wing, rho, c_wing, ...
                 Cl_alpha_2412, Cl0_2412, alpha, iw, t_wing, ...
                 "wing", ct, cr, Cm14_2412, Xc_wing);
end

% Define twist angle labels
theta_labels = arrayfun(@(x) sprintf('$\\theta=%d^\\circ$', x), theta_deg, 'UniformOutput', false);
axisFontSize = 16;
legendFontSize = 14;
lineWidth = 2;

% PLOT 1: Sectional Cl along wingspan
figure;
plot(Xc_wing(2,:), Cl_w_iso_dis', 'LineWidth', lineWidth);
xlabel('$\frac{2y}{b_w}$', 'Interpreter', 'latex', 'FontSize', axisFontSize);
ylabel('$C_l$', 'Interpreter', 'latex', 'FontSize', axisFontSize);
title('Section lift coefficient for different $\theta_t$ angles', 'Interpreter', 'latex', 'FontSize', axisFontSize);
legend(theta_labels, 'Interpreter', 'latex', 'FontSize', legendFontSize);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', axisFontSize);

% PLOT 2: Viscous Drag
figure;
plot(Xc_wing(2,:), Cd_visc_w', 'LineWidth', lineWidth);
xlabel('Wingspan [m]', 'Interpreter', 'latex', 'FontSize', axisFontSize);
ylabel('Wing $C_{D,visc}$', 'Interpreter', 'latex', 'FontSize', axisFontSize);
title('Wing $C_{D,visc}$ along wingspan', 'Interpreter', 'latex', 'FontSize', axisFontSize);
legend(theta_labels, 'Interpreter', 'latex', 'FontSize', legendFontSize);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', axisFontSize);

% PLOT 3: Induced Drag
figure;
plot(Xc_wing(2,:), Cd_induced_w', 'LineWidth', lineWidth);
xlabel('Wingspan [m]', 'Interpreter', 'latex', 'FontSize', axisFontSize);
ylabel('Wing $C_{D,ind}$', 'Interpreter', 'latex', 'FontSize', axisFontSize);
title('Wing $C_{D,ind}$ along wingspan', 'Interpreter', 'latex', 'FontSize', axisFontSize);
legend(theta_labels, 'Interpreter', 'latex', 'FontSize', legendFontSize);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', axisFontSize);

% PLOT 4: Induced AoA
figure;
plot(Xc_wing(2,:), rad2deg(alpha_ind_w)', 'LineWidth', lineWidth);
xlabel('Wingspan [m]', 'Interpreter', 'latex', 'FontSize', axisFontSize);
ylabel('Wing $\alpha_{ind}$ [$^\circ$]', 'Interpreter', 'latex', 'FontSize', axisFontSize);
title('Wing induced $\alpha$ along wingspan', 'Interpreter', 'latex', 'FontSize', axisFontSize);
legend(theta_labels, 'Interpreter', 'latex', 'FontSize', legendFontSize);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', axisFontSize);

% PLOT 5: Wing efficiency vs twist angle
figure;
plot(theta_deg, E_t_w, 'LineWidth', lineWidth);
xlabel('$\theta$ [$^\circ$]', 'Interpreter', 'latex', 'FontSize', axisFontSize);
ylabel('$E_{wing}$', 'Interpreter', 'latex', 'FontSize', axisFontSize);
title('Wing efficiency for different twist angles', 'Interpreter', 'latex', 'FontSize', axisFontSize);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', axisFontSize);

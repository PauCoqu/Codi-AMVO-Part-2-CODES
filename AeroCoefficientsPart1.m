function [Cl_alpha_0010, Cl0_0010, Cm14_0010, ...
          Cl_alpha_2412, Cl0_2412, Cm14_2412, ...
          Cl_alpha_Fowler15, Cl0_Fowler15, Cm14_Fowler15] = AeroCoefficientsPart1()

    % Angles of attack in degrees
    alpha_deg = [0, 2, 4, 6, 8];
    alpha = deg2rad(alpha_deg);

    % Lift coefficients for NACA 0010
    Cl_0010_values = [0, 0.2355, 0.4707, 0.7053, 0.9391, 1.1717];
    % Moment coefficients for NACA 0010 (about quarter-chord)
    Cm14_0010_values = [0, -0.0030, -0.0058, -0.0081, -0.0097, -0.0103];

    % Lift coefficients for NACA 2412
    Cl_2412_values = [0.2597, 0.4990, 0.7377, 0.9756, 1.2122, 1.4474];
    % Moment coefficients for NACA 2412
    Cm14_2412_values = [-0.0556, -0.0591, -0.0623, -0.0650, -0.0669, -0.0679];

    % Lift coefficients for NACA 2412 with Fowler flaps deflected 15°
    Cl_Fowler15_values = [1.4744, 3.2188, 5.2312, 0.7378, 10.0496, 12.8498];
    % Moment coefficients for NACA 2412 with Fowler flaps
    Cm14_Fowler15_values = [-0.2335, -0.2413, -0.2486, -0.2554, -0.2615, -0.2671];

    % Mean pitching moment coefficients
    Cm14_0010 = mean(Cm14_0010_values);
    Cm14_2412 = mean(Cm14_2412_values);
    Cm14_Fowler15 = mean(Cm14_Fowler15_values);

    % Linear fit (slope and intercept) for Cl vs alpha [rad] for first 5 AoA values
    p_0010 = polyfit(alpha(1:5), Cl_0010_values(1:5), 1);
    p_2412 = polyfit(alpha(1:5), Cl_2412_values(1:5), 1);
    p_Fowler15 = polyfit(alpha(1:5), Cl_Fowler15_values(1:5), 1);

    % Slopes and intercepts from the fits
    Cl_alpha_0010 = p_0010(1);
    Cl0_0010 = p_0010(2);  % Puede dejarse como 0 si se desea

    Cl_alpha_2412 = p_2412(1);
    Cl0_2412 = p_2412(2);

    Cl_alpha_Fowler15 = p_Fowler15(1);
    Cl0_Fowler15 = p_Fowler15(2);
end

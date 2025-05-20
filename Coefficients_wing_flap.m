function [Cl_t, Cm0, M_0] = ...
    Coefficients_wing_flap(ct, cr, Gamma, X, N, Xc, Q_inf, S, rho, Cm14_flap, Cm14_no_flap)

    % Inicialización de variables acumuladoras
    sum_p_Gamma = 0;
    sum_cm_coeff_1 = 0;
    sum_cm_coeff_2 = 0;
    sum_M_coeff_1 = 0;
    sum_M_coeff_2 = 0;

    % Cálculo de la cuerda media
    lambda = ct / cr;
    mean_c = (2 / 3) * cr * (1 + lambda + lambda^2) / (1 + lambda);

    % Bucle sobre los paneles
    for i = 1:N
        delta_y = abs(X(2, i + 1) - X(2, i));
        sum_p_Gamma = sum_p_Gamma + Gamma(i) * delta_y;

        if i >= round(N * 0.2) && i <= round(N * 0.8)
            sum_cm_coeff_1 = sum_cm_coeff_1 + (Xc(i) * Gamma(i) * delta_y) / (norm(Q_inf) * S * mean_c);
            sum_M_coeff_1 = sum_M_coeff_1 + (Xc(i) * Gamma(i) * delta_y);
        else
            sum_cm_coeff_2 = sum_cm_coeff_2 + (Xc(i) * Gamma(i) * delta_y) / (norm(Q_inf) * S * mean_c);
            sum_M_coeff_2 = sum_M_coeff_2 + (Xc(i) * Gamma(i) * delta_y);
        end
    end

    % Cálculo del coeficiente de sustentación total
    L = sum_p_Gamma * rho * Q_inf;
    Cl_t = L / (0.5 * rho * norm(Q_inf)^2 * S);

    % Cálculo de los coeficientes de momento
    Cm0_1 = Cm14_flap - 2 * sum_cm_coeff_1;
    Cm0_2 = Cm14_no_flap - 2 * sum_cm_coeff_2;
    Cm0 = Cm0_1 + Cm0_2;

    % Cálculo del momento total
    M_01 = 0.5 * rho * norm(Q_inf)^2 * S * mean_c * Cm14_flap - rho * norm(Q_inf) * sum_M_coeff_1;
    M_02 = 0.5 * rho * norm(Q_inf)^2 * S * mean_c * Cm14_no_flap - rho * norm(Q_inf) * sum_M_coeff_2;
    M_0 = M_01 + M_02;

end

function [Cl_t, Cl_dis, alpha_ind, Cd_visc, Cd_visc_t, Cd_induced, Cd_induced_t, ...
          E_t, Cd_t, L, Cm0, M_0] = ...
    Lift_and_CL(N, Gamma, X, Q_inf, S, rho, c, Cl_alpha, Cl_0, alpha, inc, twist, opt, ct, cr, Cm14, Xc)

    % Inicialización de variables
    sum_p_Gamma = 0;
    D_induced_t = 0;
    Cl_dis = zeros(1, N);
    Cd_visc = zeros(1, N);
    Cd_induced = zeros(1, N);
    alpha_ind = zeros(1, N);
    sum_Cd_visc = 0;

    % Bucle sobre los segmentos
    for i = 1:N
        delta_y = abs(X(2, i+1) - X(2, i));
        sum_p_Gamma = sum_p_Gamma + Gamma(i) * delta_y;

        % Coeficiente de sustentación local
        Cl_dis(i) = (2 * Gamma(i)) / (c(i) * norm(Q_inf));

        % Ángulo de ataque inducido
        alpha_ind(i) = (Cl_dis(i) - Cl_0) / Cl_alpha - (alpha + inc) - twist(i);

        % Cálculo del coeficiente de arrastre viscoso
        if opt == "wing"
            Cd_visc(i) = 0.0063 - 0.0033 * Cl_dis(i) + 0.0067 * Cl_dis(i)^2;
        elseif opt == "canard"
            Cd_visc(i) = 0.0075 + 0.0055 * Cl_dis(i)^2;
        end

        % Arrastre inducido local
        D_induced = -rho * norm(Q_inf) * Gamma(i) * delta_y * alpha_ind(i);
        Cd_induced(i) = D_induced / (0.5 * rho * norm(Q_inf)^2 * S);

        % Suma para Cd_visc total
        sum_Cd_visc = sum_Cd_visc + Cd_visc(i) * c(i) * delta_y;

        % Suma para D_induced total
        D_induced_t = D_induced_t + D_induced;
    end

    % Cálculo de coeficientes globales
    Cd_induced_t = (2 * D_induced_t) / (rho * norm(Q_inf)^2 * S);
    Cd_visc_t = sum_Cd_visc / S;
    L = sum_p_Gamma * rho * Q_inf;
    Cl_t = L / (0.5 * rho * norm(Q_inf)^2 * S);
    Cd_t = Cd_visc_t + Cd_induced_t;
    E_t = Cl_t / Cd_t;

    % Cálculo de Cm0 y M0
    lambda = ct / cr;
    mean_c = (2/3) * cr * (1 + lambda + lambda^2) / (1 + lambda);
    sum_cm_coeff = 0;
    sum_M_coeff = 0;

    for j = 1:N
        delta_y = abs(X(2, j+1) - X(2, j));
        sum_cm_coeff = sum_cm_coeff + (Xc(j) * Gamma(j) * delta_y) / (norm(Q_inf) * S * mean_c);
        sum_M_coeff = sum_M_coeff + (Xc(j) * Gamma(j) * delta_y);
    end

    Cm0 = Cm14 - 2 * sum_cm_coeff;
    M_0 = 0.5 * rho * norm(Q_inf)^2 * S * mean_c * Cm14 - rho * norm(Q_inf) * sum_M_coeff;
end

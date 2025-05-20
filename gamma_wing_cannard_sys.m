function [Gamma_w, Gamma_c, A, b] = ...
    gamma_wing_cannard_sys(Q_inf, alpha, Ur, Nw, Nc, X_wing, Xc_wing, c_wing, i_w, t_wing, ...
    Cl_alpha_2412, Cl0_2412, X_canard, Xc_canard, c_canard, i_c, Cl_alpha_0010, Cl0_0010)

    % Número total de incógnitas
    N = Nw + Nc;
    A = zeros(N, N);
    b = zeros(N, 1);
    k = [0; 0; 1]; % Dirección z para proyección

    %% WING-WING
    for i = 1:Nw
        b_ww(i,1) = 0.5 * c_wing(i) * Q_inf * (Cl0_2412 + Cl_alpha_2412 * (alpha + t_wing(i) + i_w));
        for j = 1:Nw
            if j == i
                V_ii = self_vortex(Xc_wing(:, i), X_wing(:, j), X_wing(:, j+1), Ur);
                a_ww(i, j) = -0.5 * Cl_alpha_2412 * c_wing(j) * dot(V_ii, k) + 1;
            else
                V_ij = horseshoe_vortex(Xc_wing(:, i), X_wing(:, j), X_wing(:, j+1), Ur);
                a_ww(i, j) = -0.5 * Cl_alpha_2412 * c_wing(i) * dot(V_ij, k);
            end
        end
    end

    %% CANARD-WING (efecto del canard en el ala)
    for i = 1:Nw
        for j = 1:Nc
            V_ij_comp = horseshoe_vortex(Xc_wing(:, i), X_canard(:, j), X_canard(:, j+1), Ur);
            a_cw(i, j) = -0.5 * Cl_alpha_2412 * c_wing(i) * dot(V_ij_comp, k);
        end
    end

    %% CANARD-CANARD
    for i = 1:Nc
        b_cc(i,1) = 0.5 * c_canard(i) * Q_inf * (Cl0_0010 + Cl_alpha_0010 * (alpha + i_c)); % sin twist en el canard
        for j = 1:Nc
            if j == i
                V_ii = self_vortex(Xc_canard(:, j), X_canard(:, j), X_canard(:, j+1), Ur);
                a_cc(i, j) = -0.5 * Cl_alpha_0010 * c_canard(j) * dot(V_ii, k) + 1;
            else
                V_ij = horseshoe_vortex(Xc_canard(:, i), X_canard(:, j), X_canard(:, j+1), Ur);
                a_cc(i, j) = -0.5 * Cl_alpha_0010 * c_canard(i) * dot(V_ij, k);
            end
        end
    end

    %% WING-CANARD (efecto del ala en el canard)
    for i = 1:Nc
        for j = 1:Nw
            V_ij_comp2 = horseshoe_vortex(Xc_canard(:, i), X_wing(:, j), X_wing(:, j+1), Ur);
            a_wc(i, j) = -0.5 * Cl_alpha_0010 * c_canard(i) * dot(V_ij_comp2, k);
        end
    end

    %% Ensamblaje del sistema
    A(1:Nw, 1:Nw) = a_ww;
    A(Nw+1:N, Nw+1:N) = a_cc;
    A(1:Nw, Nw+1:N) = a_cw;
    A(Nw+1:N, 1:Nw) = a_wc;

    b(1:Nw) = b_ww;
    b(Nw+1:N) = b_cc;

    %% Resolución del sistema
    Gamma = A \ b;

    % Separación de circulaciones
    Gamma_w = Gamma(1:Nw);
    Gamma_c = Gamma(Nw+1:N);

end

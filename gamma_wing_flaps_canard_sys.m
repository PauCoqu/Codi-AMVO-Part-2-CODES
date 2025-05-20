function [Gamma_w_4, Gamma_c_4] = ...
    gamma_wing_flaps_canard_sys(Q_inf, alpha, Ur, Nw, Nc, ...
    X_wing, Xc_wing, c_wing, i_w, t_wing, ...
    Cl_alpha_2412, Cl0_2412, ...
    X_canard, Xc_canard, c_canard, i_c, ...
    Cl_alpha_0010, Cl0_0010, ...
    Cl0_Fowler15, Cl_alpha_Fowler15)

    % Dimensiones del sistema
    N = Nw + Nc;
    A = zeros(N, N);
    b = zeros(N, 1);
    k = [0; 0; 1];

    %% WING-WING
    for i = 1:Nw
        if i >= round(Nw * 0.2) && i <= round(Nw * 0.8)
            b_ww(i,1) = 0.5 * c_wing(i) * Q_inf * (Cl0_Fowler15 + Cl_alpha_Fowler15 * (alpha + t_wing(i) + i_w));
            for j = 1:Nw
                if j == i
                    V_ii = self_vortex(Xc_wing(:, i), X_wing(:, j), X_wing(:, j + 1), Ur);
                    a_ww(i, j) = -0.5 * Cl_alpha_Fowler15 * c_wing(i) * dot(V_ii, k) + 1;
                else
                    V_ij = horseshoe_vortex(Xc_wing(:, i), X_wing(:, j), X_wing(:, j + 1), Ur);
                    a_ww(i, j) = -0.5 * Cl_alpha_Fowler15 * c_wing(i) * dot(V_ij, k);
                end
            end
        else
            b_ww(i,1) = 0.5 * c_wing(i) * Q_inf * (Cl0_2412 + Cl_alpha_2412 * (alpha + t_wing(i) + i_w));
            for j = 1:Nw
                if j == i
                    V_ii = self_vortex(Xc_wing(:, i), X_wing(:, j), X_wing(:, j + 1), Ur);
                    a_ww(i, j) = -0.5 * Cl_alpha_2412 * c_wing(j) * dot(V_ii, k) + 1;
                else
                    V_ij = horseshoe_vortex(Xc_wing(:, i), X_wing(:, j), X_wing(:, j + 1), Ur);
                    a_ww(i, j) = -0.5 * Cl_alpha_2412 * c_wing(i) * dot(V_ij, k);
                end
            end
        end
    end

    %% CANARD-WING
    for i = 1:Nw
        if i >= round(Nw * 0.2) && i <= round(Nw * 0.8)
            for j = 1:Nc
                V_ij_comp = horseshoe_vortex(Xc_wing(:, i), X_canard(:, j), X_canard(:, j + 1), Ur);
                a_cw(i, j) = -0.5 * Cl_alpha_Fowler15 * c_wing(i) * dot(V_ij_comp, k);
            end
        else
            for j = 1:Nc
                V_ij_comp = horseshoe_vortex(Xc_wing(:, i), X_canard(:, j), X_canard(:, j + 1), Ur);
                a_cw(i, j) = -0.5 * Cl_alpha_2412 * c_wing(i) * dot(V_ij_comp, k);
            end
        end
    end

    %% CANARD-CANARD
    for i = 1:Nc
        b_cc(i,1) = 0.5 * c_canard(i) * Q_inf * (Cl0_0010 + Cl_alpha_0010 * (alpha + i_c)); % sin twist
        for j = 1:Nc
            if j == i
                V_ii = self_vortex(Xc_canard(:, j), X_canard(:, j), X_canard(:, j + 1), Ur);
                a_cc(i, j) = -0.5 * Cl_alpha_0010 * c_canard(j) * dot(V_ii, k) + 1;
            else
                V_ij = horseshoe_vortex(Xc_canard(:, i), X_canard(:, j), X_canard(:, j + 1), Ur);
                a_cc(i, j) = -0.5 * Cl_alpha_0010 * c_canard(i) * dot(V_ij, k);
            end
        end
    end

    %% WING-CANARD
    for i = 1:Nc
        for j = 1:Nw
            V_ij_comp2 = horseshoe_vortex(Xc_canard(:, i), X_wing(:, j), X_wing(:, j + 1), Ur);
            a_wc(i, j) = -0.5 * Cl_alpha_0010 * c_canard(i) * dot(V_ij_comp2, k);
        end
    end

    %% Ensamblado del sistema
    A(1:Nw, 1:Nw) = a_ww;
    A(Nw+1:N, Nw+1:N) = a_cc;
    A(1:Nw, Nw+1:N) = a_cw;
    A(Nw+1:N, 1:Nw) = a_wc;

    b(1:Nw) = b_ww;
    b(Nw+1:N) = b_cc;

    %% Solución del sistema
    Gamma = A \ b;

    % Separación
    Gamma_w_4 = Gamma(1:Nw);
    Gamma_c_4 = Gamma(Nw+1:N);
end

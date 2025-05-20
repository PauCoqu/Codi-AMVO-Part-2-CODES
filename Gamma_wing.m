function [gamma_isol] = Gamma_wing(N, Q_inf, alpha, Ur, X, Xc, c, inc, t, Cl_alpha, Cl0)
    % Inicialización de matrices del sistema lineal
    a_ww = zeros(N, N);
    b_ww = zeros(N, 1);

    % Vector unitario perpendicular al plano de la corriente
    k = [-sin(alpha + inc); 0; cos(alpha + inc)];

    % Construcción de la matriz A y del vector b
    for i = 1:N
        % Término constante para el sistema
        b_ww(i) = 0.5 * c(i) * Q_inf * (Cl0 + Cl_alpha * (alpha + t(i) + inc));

        for j = 1:N
            if j == i
                % Velocidad inducida por el propio vórtice (self-induced)
                V_ii = self_vortex(Xc(:, i), X(:, j), X(:, j + 1), Ur);
                a_ww(i, j) = -0.5 * Cl_alpha * c(i) * dot(V_ii, k) + 1;
            else
                % Velocidad inducida por otros vórtices (interacción)
                V_ij = horseshoe_vortex(Xc(:, i), X(:, j), X(:, j + 1), Ur);
                a_ww(i, j) = -0.5 * Cl_alpha * c(i) * dot(V_ij, k);
            end
        end
    end

    % Resolución del sistema lineal A * gamma = b
    gamma_isol = a_ww \ b_ww;
end

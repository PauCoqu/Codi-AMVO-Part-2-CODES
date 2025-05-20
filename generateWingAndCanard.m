function [X_wing, Xc_wing, c_wing, X_canard, Xc_canard, c_canard] = ...
    generateWingAndCanard(Nw, b, cr, ct, Nc, bh, crh, cth, lh, theta)

    % ========================
    % Wing Geometry
    % ========================

    % Coordenadas del ala (alineada con el eje x)
    wing_x = zeros(1, Nw + 1);                % eje x (longitudinal)
    wing_y = linspace(-b/2, b/2, Nw + 1);     % eje y (envergadura)
    wing_z = zeros(1, Nw + 1);                % eje z (altura)
    X_wing = [wing_x; wing_y; wing_z];        % matriz de posiciones de los vértices

    % Puntos de control en el centro de cada segmento
    Xc_wing = (X_wing(:, 1:Nw) + X_wing(:, 2:Nw + 1)) / 2;

    % Distribución de la cuerda del ala
    c_wing = ComputeChordDistribution(cr, ct, b, Xc_wing);

    % ========================
    % Canard Geometry
    % ========================

    % Coordenadas del canard (desplazado en x, más alto en z)
    canard_x = -lh * ones(1, Nc + 1);                 % eje x (adelantado respecto al ala)
    canard_y = linspace(-bh/2, bh/2, Nc + 1);         % eje y (envergadura)
    canard_z = 0.1 * ones(1, Nc + 1);                 % eje z (ligeramente elevado)
    X_canard = [canard_x; canard_y; canard_z];        % matriz de posiciones del canard

    % Puntos de control en el centro de cada segmento del canard
    Xc_canard = (X_canard(:, 1:Nc) + X_canard(:, 2:Nc + 1)) / 2;

    % Distribución de la cuerda del canard
    c_canard = ComputeChordDistribution(crh, cth, bh, Xc_canard);
end

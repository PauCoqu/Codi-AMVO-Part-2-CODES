function [V_ij] = horseshoe_vortex(X_c, X_1, X_2, Ur)
% Calculates the induced velocity at control point X_c due to a horseshoe vortex
% formed by bound vortex segment from X_1 to X_2 and two trailing semi-infinite vortices.

    % Vectores desde el punto de control hacia los extremos del segmento de vórtice
    r1 = X_c - X_1;
    r2 = X_c - X_2;

    % Velocidades inducidas por los extremos infinitos del vórtice
    V_infA = SemiInfVortex(r1, Ur);
    V_infB = SemiInfVortex(r2, Ur);

    % Velocidad inducida por el segmento finito superior (entre X_1 y X_2)
    V_AB = TopVortex(r1, r2);

    % Resultado total: suma de las contribuciones
    V_ij = V_infA + V_AB - V_infB;
end

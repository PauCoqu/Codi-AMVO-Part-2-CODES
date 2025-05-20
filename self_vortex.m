function [V_ii] = self_vortex(X_c, X_1, X_2, Ur)
% SELF_VORTEX Calcula la velocidad inducida por un vórtice en sí mismo
% Entrada:
%   X_c  - Punto de control (3x1)
%   X_1  - Extremo delantero del vórtice (3x1)
%   X_2  - Extremo trasero del vórtice (3x1)
%   Ur   - Dirección del flujo libre (3x1)
% Salida:
%   V_ii - Velocidad inducida en el punto X_c por el propio vórtice

    r1 = X_c - X_1;
    r2 = X_c - X_2;

    V_infA = SemiInfVortex(r1, Ur);
    V_infB = SemiInfVortex(r2, Ur);

    V_ii = V_infA - V_infB;
end

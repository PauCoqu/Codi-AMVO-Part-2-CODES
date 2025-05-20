function [t_wing] = ComputeTwistDistribution(theta_w, Nw)
% COMPUTETWISTDISTRIBUTION Genera la distribución de twist simétrica
% Entrada:
%   theta_w - Ángulo de twist máximo (en radianes)
%   Nw - Número de segmentos en la envergadura
% Salida:
%   t_wing - Vector de twist en cada segmento

    % Distribución simétrica de twist en la semienvergadura
    t_wing_1 = linspace(theta_w, 0, Nw/2 + 1);       % desde raíz hasta centro
    t_wing_2 = linspace(0, theta_w, Nw/2 + 1);       % desde centro hasta la punta

    % Eliminar el punto duplicado del centro (ya incluido en ambos)
    t_wing_n = [t_wing_1, t_wing_2(2:end)];

    % Promediar en cada panel (valor medio entre dos puntos consecutivos)
    for i = 1:Nw
        t_wing(i) = (t_wing_n(i) + t_wing_n(i + 1)) / 2;
    end

    % Convertir a vector fila (opcional, dependiendo de uso posterior)
    t_wing = t_wing(:)';
end

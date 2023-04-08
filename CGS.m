% Parámetros
N = 100;
a = 0.047;
o_values = linspace(pi/2, 3*pi/2, N+1);
o_values(end) = [];

% Calcular L y Cg
[Cg_values, L_values] = calc_L_Cg(o_values, a);

% Funciones de interpolación    
L_spline = spline(o_values, L_values);
Cg_spline = spline(o_values, Cg_values);

% Evaluar splines
o_eval = linspace(pi/2, 3*pi/2, 1000);
L_eval = ppval(L_spline, o_eval);
Cg_eval = ppval(Cg_spline, o_eval);

% Gráficas
figure;
plot(o_eval, L_eval);
title('L en función de \theta');
xlabel('\theta');
ylabel('L');

figure;
plot(o_eval, Cg_eval);
title('Cg en función de \theta');
xlabel('\theta');
ylabel('Cg');

% Función para calcular L y Cg
function [Cg_values, L_values] = calc_L_Cg(o_values, a)
    L_values = zeros(size(o_values)); % inicialización de valores de L (longitud de arco)
    Cg_values = zeros(size(o_values)); % inicialización de valores de Cg

    for i = 1:length(o_values)
        o = o_values(i);
        t_values = linspace(pi/2, 3*pi/2, 1000);
        x = a .* (((sin(o) .* cos(t_values - o) + sin(t_values)) ./ (1 + sin(t_values - o))) + cos(o));
        y = -a .* (((cos(o) .* cos(t_values - o) + sin(t_values)) ./ (1 + sin(t_values - o))) - sin(o));
        dx = diff(x);
        dy = diff(y);
        L = sum(sqrt(dx.^2 + dy.^2));
        Cg = (2 .* x(end)) / (2 .* pi .* a);

        L_values(i) = L;
        Cg_values(i) = Cg;
    end
end

% Nuevos valores de ángulo y colores
angulos = [pi/2, pi, 3*pi/2];
colores = {'red', 'green', 'blue'};

% Gráfica 1: Sv vs Cg
o_values = linspace(pi/2, 3*pi/2, 100);
a = 0.047;
[Cg_values, Sv_values] = calc_Sv_Cg(o_values, a);

figure(1)
plot(Cg_values, Sv_values)
xlabel('Cg')
ylabel('Sv')
title('Esbeltez vs CG')

% Líneas verticales en Gráfica 1: Sv vs Cg
for i = 1:length(angulos)
    plot_red_line(angulos(i), colores{i}, a);
end

% Función para calcular Sv y Cg
function [Cg_values, Sv_values] = calc_Sv_Cg(o_values, a)
    Sv_values = zeros(size(o_values));
    Cg_values = zeros(size(o_values));

    for i = 1:length(o_values)
        o = o_values(i);
        t = linspace(pi/2, 3*pi/2, 100);
        n = length(t);

        x = a * (((sin(o) * cos(t - o) + n * cos(t)) ./ (1 + sin(t - o))) + cos(o));
        y = -a * (((cos(o) * cos(t - o) + n * cos(t)) ./ (1 + sin(t - o))) - sin(o));

        Sv = y ./ (2 * x - a);
        Cg = (2 * x) ./ (2 * pi * a);

        Sv_values(i) = mean(Sv);
        Cg_values(i) = mean(Cg);
    end
end

% Función para trazar las líneas verticales en la gráfica 1
function plot_red_line(o_fixed, color, a)
    t_values = linspace(pi/2, 3*pi/2, 1000);
    n = length(t_values);

    Sv_values = zeros(size(t_values));
    Cg_values = zeros(size(t_values));

    for i = 1:length(t_values)
        t = t_values(i);
        x = a * (((sin(o_fixed) * cos(t - o_fixed) + n * cos(t)) ./ (1 + sin(t - o_fixed))) + cos(o_fixed));
        y = -a * (((cos(o_fixed) * cos(t - o_fixed) + n * cos(t)) ./ (1 + sin(t - o_fixed))) - sin(o_fixed));

        Sv = y ./ (2 * x - a);
        Cg = (2 * x) ./ (2 * pi * a);

        Sv_values(i) = Sv;
        Cg_values(i) = Cg;
    end

    hold on
    plot(Cg_values, Sv_values, color)
    hold off
end

o_values = linspace(0, pi/2, 40); % valores de "o" de 0 a pi/2
angulos = [pi/6, pi/12, pi/18, 7*pi/180, pi/36, pi/60];
colores = {'red', 'green', 'blue', 'black', 'magenta', 'yellow'};

a = 0.047;

% Gráfica 1: Sv vs Cg
[Cg_values, Sv_values] = calc_Sv_Cg(o_values, a);
figure(1)
plot(Cg_values, Sv_values)
xlabel('Cg')
ylabel('Sv')
title('Esbeltez vs CG')

for i = 1:length(angulos)
    plot_red_line(angulos(i), colores{i}, a);
end

% Líneas verticales en Gráfica 1: Sv vs Cg
figure(1)
for i = 1:length(angulos)
    o = angulos(i);
    t = linspace(pi/2, 3*pi/2, 100);
    x = a*(((sin(o).*cos(t-o)+1.*cos(t))./(1+sin(t-o)))+cos(o));
    y = -a*(((cos(o).*cos(t-o)+1.*cos(t))./(1+sin(t-o)))-sin(o));
    Cg = (2*x)./(2*pi*a);
    plot_vertical_line(Cg, colores{i});
end

% Gráfica 2: Longitud de arco vs Cg
[Cg_values, L_values] = calc_L_Cg(o_values, a);
figure(2)
plot(Cg_values, L_values)
xlabel('Cg')
ylabel('L (longitud de arco)')
title('Longitud de arco vs Cg')

% Líneas verticales en Gráfica 2: Longitud de arco vs Cg
figure(2)
for i = 1:length(angulos)
    o = angulos(i);
    t = linspace(pi/2, 3*pi/2, 100);
    x = a*(((sin(o).*cos(t-o)+1.*cos(t))./(1+sin(t-o)))+cos(o));
    y = -a*(((cos(o).*cos(t-o)+1.*cos(t))./(1+sin(t-o)))-sin(o));
    Cg = (2*x)./(2*pi*a);
    plot_vertical_line(Cg, colores{i});
end


% Cambiar el valor de "a" en las funciones correspondientes
a_value = 0.047;
[Cg_values, Sv_values] = calc_Sv_Cg(o_values, a_value);
[Cg_values, L_values] = calc_L_Cg(o_values, a_value);

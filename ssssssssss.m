o_values = linspace(0, pi/2, 40); % valores de "o" de 0 a pi/2
angulos = [pi/6, pi/12, pi/18, 7*pi/180, pi/36, pi/60];
colores = {'red', 'green', 'blue', 'black', 'magenta', 'yellow'};

% Gráfica 1: Sv vs Cg
[Cg_values, Sv_values] = calc_Sv_Cg(o_values);
figure(1)
plot(Cg_values, Sv_values)
xlabel('Cg')
ylabel('Sv')
title('Esbeltez vs CG')

for i = 1:length(angulos)
    plot_red_line(angulos(i), colores{i});
end

% Líneas verticales en Gráfica 1: Sv vs Cg
figure(1)
xline(((2*(0.047*(1+sin(pi/6)).*cos(pi/2 - pi/6))./(1-sin(pi/2 - 2*pi/6)))-0.047)/0.047, 'red','30º'); % o = pi/6
xline(((2*(0.047*(1+sin(pi/12)).*cos(pi/2 - pi/12))./(1-sin(pi/2 - 2*pi/12)))-0.047)/0.047, 'green','15º');% o = pi/12
xline(((2*(0.047*(1+sin(pi/18)).*cos(pi/2 - pi/18))./(1-sin(pi/2 - 2*pi/18)))-0.047)/0.047, 'blue','10º'); % o = pi/18
xline(((2*(0.047*(1+sin(pi*7/180)).*cos(pi/2 - pi*7/180))./(1-sin(pi/2 - 2*pi*7/180)))-0.047)/0.047, 'black','7º'); % o = pi*7/180
xline(((2*(0.047*(1+sin(pi/36)).*cos(pi/2 - pi/36))./(1-sin(pi/2 - 2*pi/36)))-0.047)/0.047, 'magenta','5º'); % o = pi/36
xline(((2*(0.047*(1+sin(pi/60)).*cos(pi/2 - pi/60))./(1-sin(pi/2 - 2*pi/60)))-0.047)/0.047, 'yellow','3º'); % o = pi/60

o_values = angulos; % valores de ángulo fijo

% Gráfica 2: Longitud de arco vs Cg
[Cg_values, L_values] = calc_L_Cg(o_values);
figure(2)
plot(Cg_values, L_values)
xlabel('Cg')
ylabel('L (longitud de arco)')
title('Longitud de arco vs Cg')

% Líneas verticales en Gráfica 2: Longitud de arco vs Cg
figure(2)
xline(((2*(0.047*(1+sin(pi/6)).*cos(pi/2 - pi/6))./(1-sin(pi/2 - 2*pi/6)))-0.047)/0.047, 'r-', '30º'); % o = pi/6
xline(((2*(0.047*(1+sin(pi/12)).*cos(pi/2 - pi/12))./(1-sin(pi/2 - 2*pi/12)))-0.047)/0.047, 'g-', '15º'); % o = pi/12
xline(((2*(0.047*(1+sin(pi/18)).*cos(pi/2 - pi/18))./(1-sin(pi/2 - 2*pi/18)))-0.047)/0.047, 'b-', '10º'); % o = pi/18
xline(((2*(0.047*(1+sin(pi*7/180)).*cos(pi/2 - pi*7/180))./(1-sin(pi/2 - 2*pi*7/180)))-0.047)/0.047, 'k-', '7º'); % o = pi*7/180
xline(((2*(0.047*(1+sin(pi/36)).*cos(pi/2 - pi/36))./(1-sin(pi/2 - 2*pi/36)))-0.047)/0.047, 'm-', '5º'); % o = pi/36
xline(((2*(0.047*(1+sin(pi/60)).*cos(pi/2 - pi/60))./(1-sin(pi/2 - 2*pi/60)))-0.047)/0.047, 'y-', '3º'); % o = pi/60


% Líneas verticales en Gráfica 2: Longitud de arco vs Cg
xline(Cg_values(1), 'k--', 'Cg = 0');
for i = 2:length(o_values)
    o_fixed = o_values(i);
    [Cg_values_fixed, L_values_fixed] = calc_L_Cg_fixed_o(o_fixed);
    %Cg_intersect = interp1(L_values_fixed, Cg_values_fixed, 0.0001);
    %xline(Cg_intersect, 'k--', ['Cg = ' num2str(Cg_intersect)]);
end


function [Cg_values, Sv_values] = calc_Sv_Cg(o_values)
a=0.5;
    Sv_values = zeros(size(o_values)); % inicialización de valores de Sv
    Cg_values = zeros(size(o_values)); % inicialización de valores de Cg

    for i = 1:length(o_values)
        o = o_values(i);
        t = pi/2 - o;
        x = a*((sin(o)*cos(t-o)+1*cos(t))/(1+sin(t-o))+cos(o));
        y = -a*((cos(o)*cos(t-o)+1*cos(t))/(1+sin(t-o))-sin(o));
        Sv = y/(2*x-a);
        Cg = (2*x)/(2*pi*a);
        Sv_values(i) = Sv;
        Cg_values(i) = Cg;
    end
end


function [Cg_values, L_values] = calc_L_Cg_fixed_o(o_values)
    L_values = zeros(size(o_values)); % inicialización de valores de L (longitud de arco)
    Cg_values = zeros(size(o_values)); % inicialización de valores de Cg

    for i = 1:length(o_values)
        o = o_values(i);
        tmax = pi/2 - o;
        
        dy_dt = @(t) (0.047*(1+sin(o))*(sin(o)*sin(t)-sin(t)))./(1-sin(t-o)).^2;
        dx_dt = @(t) (0.047*(1+sin(o))*(sin(o)*cos(t)-cos(t)))./(1-sin(t-o)).^2;
        
        t_values = linspace(0, tmax, 100);
        L = sum(sqrt(arrayfun(dy_dt, t_values).^2 + arrayfun(dx_dt, t_values).^2)) * (tmax / 99);
        
        Cg = (2*sin(o) - 0.5) / (0.5*sin(o));
        fprintf('o=%f, Cg=%f, L=%f\n', o, Cg, L);

        L_values(i) = L;
        Cg_values(i) = Cg;
        
        % Imprimir los valores de Cg y L
        fprintf('o=%f, Cg=%f, L=%f\n', o, Cg, L);
    end
end

% Función para calcular L y Cg
function [Cg_values, L_values] = calc_L_Cg(o_values)
    L_values = zeros(size(o_values)); % inicialización de valores de L (longitud de arco)
    Cg_values = zeros(size(o_values)); % inicialización de valores de Cg

    for i = 1:length(o_values)
        o = o_values(i);
        tmax = pi/2 - o;
        
        %y_fun = @(t) (0.047*(1+sin(o))*sin(t))/(1-sin(t-o));
        %x_fun = @(t) (0.047*(1+sin(o))*cos(t))/(1-sin(t-o));
        
        dy_dt = @(t) (0.047*(1+sin(o))*(sin(o)*sin(t)-sin(t)))./(1-sin(t-o)).^2;
        dx_dt = @(t) (0.047*(1+sin(o))*(sin(o)*cos(t)-cos(t)))./(1-sin(t-o)).^2;
        
        t_values = linspace(0, tmax, 100);
        L = sum(sqrt(arrayfun(dy_dt, t_values).^2 + arrayfun(dx_dt, t_values).^2)) * (tmax / 99);
        
        Cg = 1 / sin(o);
        
        L_values(i) = L;
        Cg_values(i) = Cg;
    end
end





% Función para calcular las líneas rojas de la gráfica 1 (Sv vs Cg)
function plot_red_line(o_fixed, color)
    tmax = pi/2 - o_fixed;
    t_values = linspace(0, tmax, 1000);
    
    Sv_values = zeros(size(t_values));
    Cg_values = zeros(size(t_values));
    
    for i = 1:length(t_values)
        t = t_values(i);
        y = (0.5*(1+sin(o_fixed))*sin(t))/(1-sin(t-o_fixed));
        x = (0.5*(1+sin(o_fixed))*cos(t))/(1-sin(t-o_fixed));
        Sv = y/(2*x-0.5);
        Cg = (2*x-0.5)/0.5;
        
                Sv_values(i) = Sv;
        Cg_values(i) = Cg;
    end
    
    hold on
    plot(Cg_values, Sv_values, color)
    hold off
end

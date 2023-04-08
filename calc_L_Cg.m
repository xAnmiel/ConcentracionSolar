function [Cg_values, L_values] = calc_L_Cg(o_values, a)
    L_values = zeros(size(o_values));
    Cg_values = zeros(size(o_values));

    for i = 1:length(o_values)
        o = o_values(i);
        tmax = 3*pi/2 - o;
        
        dy_dt = @(t) (-a*(cos(o)*cos(t-o)-sin(o))*(sin(t-o)+1))./((sin(t-o)+1).^2) - a*sin(o);
        dx_dt = @(t) (a*(sin(o)*cos(t-o)+cos(t))*(sin(t-o)+1))./((sin(t-o)+1).^2) + a*cos(o);
        
        t_values = linspace(pi/2, tmax, 100);
        L = sum(sqrt(arrayfun(dy_dt, t_values).^2 + arrayfun(dx_dt, t_values).^2)) * ((tmax - pi/2) / 99);
        
        Cg = (2*sin(o) - 0.047) / (0.047*sin(o));
        
        L_values(i) = L;
        Cg_values(i) = Cg;
    end
end

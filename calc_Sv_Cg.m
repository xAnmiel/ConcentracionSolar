function [Cg_values, Sv_values] = calc_Sv_Cg(o_values, a)
    Sv_values = zeros(size(o_values));
    Cg_values = zeros(size(o_values));

    for i = 1:length(o_values)
        o = o_values(i);
        t = 3*pi/2 - o;
        x = a*(((sin(o)*cos(t-o)+1*cos(t))/(1+sin(t-o)))+cos(o));
        y = -a*(((cos(o)*cos(t-o)+1*cos(t))/(1+sin(t-o)))-sin(o));
        Sv = y/(2*x);
        Cg = (2*x)/(2*pi*a);
        Sv_values(i) = Sv;
        Cg_values(i) = Cg;
    end
end

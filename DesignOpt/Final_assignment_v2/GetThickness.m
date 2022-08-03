function [t_c] = GetThickness(xc)
naca = '0008';
t = str2double(naca(3:4)) / 100;

y_c = 5 * t * (0.2969 * sqrt(xc) - 0.1260 * xc - ...
    0.3516 * xc^2 + 0.2843 * xc^3 - 0.1036 * xc^4);

t_c = 2 * y_c;

end



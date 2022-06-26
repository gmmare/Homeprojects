x = [0.1000    0.6000    0.0015    0.0015    0.0008    0.0009];
hi = 1e-8;
dg1dx = Sensitivity_analysis(xq, hi)

function [dg1dx] = Sensitivity_analysis(xq, hi)
opt_params;
%setting change in x
dI_xx = zeros(1, 6);

%flutter value 
Ixx_local = GetInertia(xq);
p_max1 = constraints(Ixx_local);
p_max2 = constraints(Ixx_local + hi);

% constraint value
g1x = 1 - p_max1/p_ref;
g1xdx = 1 - p_max2/p_ref;

dg1dIxx = (g1xdx - g1x)/hi;


for i=1:length(xq)

    %setting up new x vector
    dx = zeros(1, 6);
    dx(i) =  hi;

    x_new = xq + dx;
    %calculating new objective value
    Ixx_new = GetInertia(x_new);
    dI_xx(i) = (Ixx_new - Ixx_local)/hi;
end

dg1dx = dg1dIxx * dI_xx;

end



function [f, b1, b2, g1, g2] = objective(x)
opt_params;

%analysis tools
W = GetWeight(x);
A = GetVolume(x);
I_xx = GetInertia(x);

%objective value
f_obj = c * (1-k)*(A_ref/A) + k*(W/W_ref);
%flutter constraint
p_max = constraints(I_xx);

%constraint forumlation
g1 = 1 - p_max/p_ref;
if p_max > p_ref
    g1 = abs(p_max/p_ref);
end

g2 = 0.8 * I_xx_ref/(I_xx) - 1;

mult_barrier = 12;
% mult_barrier = 0;
%calculating barrier function
b1 = real( - (1/r) * log(-g1));
b2 = mult_barrier *real( - (1/r) * log(-g2));

f = f_obj;

end


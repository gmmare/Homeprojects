function [f, b1, b2] = objective(x)
opt_params;

%analysis tools
W = GetWeight(x);
A = GetVolume(x);

%objective value
f_obj = (1-k)*(A/A_ref) + k*(W/W_ref);

%flutter constraint
p_max = constraints(x);

%heave constraint
I_xx = GetInertia(x);


%constraint forumlation
g1 = 1 - p_max/p_ref;
if p_max > p_ref
    g1 = abs(p_max/p_ref);
end

g2 = I_xx/I_xx_ref - 1;


%calculating barrier function
b1 = real( - (1/r) * log(-g1));
b2 = real( - (1/r) * log(-g2));

f = f_obj;

end


% xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];
% objectivet(xq)
% % 

function [f, b1, b2, b3] = objective(x)
opt_params;

%analysis tools
W = GetWeight(x);
A = GetVolume(x);

%objective value
f_obj = (1-k)*(A/A_ref) + k*(W/W_ref);

%flutter constraint
p_max = constraints(x);

%heave constraint
h_max = constraint2(x);

%constraint forumlation
g1 = 1 - p_max/p_ref;
g2 = h_max/h_ref - 1;
g3 = 0.6 * (h_ref/h_max) - 1;

%calculating barrier function
b1 = real(-(1/r) * log(-g1));
b2 = real( - (1/r) * log(-g2));
b3 = real( - (1/r) * log(-g3));


f = f_obj;

end


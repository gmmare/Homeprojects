function [W] = GetWeight(x)
%GETWEIGHT calculate the weight of the structural components. Make
%simplification of an I beam
opt_params;

%{ 
design vector x
1 Front spar location (continuous)
2 Aft spar location (continuous)
3 Web-thickness front spar (continuous)
4 Web-thickness aft spar (continuous)
5 Flange thickness front spar (continuous)
6 Flange thickness aft spar (continuous)
%}

%Front spar
w_web_f = GetThickness(x(1)) * x(3);
w_flange_f = 2 * x(5) * flange_width;

%aft spar
w_web_a = GetThickness(x(2)) * x(4);
w_flange_a = 2 * x(6) * flange_width;

W = w_web_f + w_web_a + w_flange_f + w_flange_a;

end


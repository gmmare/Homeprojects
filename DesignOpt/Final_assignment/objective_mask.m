function [f_sum] = objective_mask(alpha, xq, sq)
%getting optimizer parameters
opt_params;

% line searching
x_des = xq + alpha*sq; 

%objective value
[f, b1, b2, b3] = objective(x_des);
f_sum = f + b1;% + b2 + b3;
end


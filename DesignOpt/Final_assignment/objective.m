function [f] = objective(x)
%getting optimizer parameters
opt_params;
x_des = zeros(1, 6);

% %de-normalize design vector first
for i= 1:length(x)
    x_des(i) = x(i) * (bounds_upper(i) - bounds_lower(i)) + bounds_lower(i);
end

f = x_des;
A_ref
I_xx_ref
end


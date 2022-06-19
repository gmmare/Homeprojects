function [I_xx] = GetInertia(x)
%INERTIA 
% design vector x
% 1 Front spar location (continuous)
% 2 Aft spar location (continuous)
% 3 Web-thickness front spar (continuous)
% 4 Web-thickness aft spar (continuous)
% 5 Flange thickness front spar (continuous)
% 6 Flange thickness aft spar (continuous)
%   Detailed explanation goes here

% flange width:
flange_width = 0.042;

% front spar
t_spar_front = GetThickness(x(1));
web_front = x(3) * t_spar_front^3 /12;
flange_front = flange_width * x(5) * (t_spar_front/2)^2;
I_xx_front = web_front + 2 * flange_front;

%aft spar
t_spar_aft = GetThickness(x(2));
web_aft = x(4) * t_spar_aft^3 /12;
flange_aft = flange_width * x(6) * (t_spar_aft/2)^2;
I_xx_aft = web_aft + 2 * flange_aft;

% total inertia
I_xx = (I_xx_front + I_xx_aft);
end


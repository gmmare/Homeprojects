% xq = [0.05 0.80 0.0012 0.0012 0.0015 0.0015];
% GetVolumet(xq)
% xq = [0.25 0.35 0.0012 0.0012 0.0015 0.0015];
% GetVolumet(xq)

function [A_total] = GetVolume(x)
% design vector x
% 1 Front spar location (continuous)
% 2 Aft spar location (continuous)
% 3 Web-thickness front spar (continuous)
% 4 Web-thickness aft spar (continuous)
% 5 Flange thickness front spar (continuous)
% 6 Flange thickness aft spar (continuous)
%   Detailed explanation goes here

% % calculate tank crosssectional area by two trapezoids
% x_mid = 0.5 * (x(1) + x(2));
% 
% % Getting thickness at each point
% h1 = GetThickness(x(1));
% h2 = GetThickness(x(2));
% hmid = GetThickness(x_mid);
% 
% % calculating trapezoidal areas
% A1 = 0.5 * abs(x_mid - x(1)) * (h1 + hmid);
% A2 = 0.5 * abs(x(2) - x_mid) * (h2 + hmid);

n_points = 40;
x_tank = linspace(x(1), x(2), n_points);
dx = abs(x(2) - x(1))/n_points;
A_total = 0;
for i=1:n_points - 1
    A_total = A_total + GetThickness(x_tank(i) + 0.5*dx) * dx;
end

end


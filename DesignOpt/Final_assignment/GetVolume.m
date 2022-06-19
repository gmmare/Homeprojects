function [A_total] = GetVolume(x)
% design vector x
% 1 Front spar location (continuous)
% 2 Aft spar location (continuous)
% 3 Web-thickness front spar (continuous)
% 4 Web-thickness aft spar (continuous)
% 5 Flange thickness front spar (continuous)
% 6 Flange thickness aft spar (continuous)
%   Detailed explanation goes here

% calculate tank crosssectional area by two trapezoids
x_mid = 0.5 * (x(1) + x(2));

% Getting thickness at each point
h1 = GetThickness(x(1));
h2 = GetThickness(x(2));
hmid = GetThickness(x_mid);

% calculating trapezoidal areas
A1 = 0.5 * abs(x_mid - x(2)) * (h1 + hmid);
A2 = 0.5 * abs(x(2) - x_mid) * (h2 + hmid);

A_total = A1 + A2;
end


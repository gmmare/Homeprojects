xq = [0.1674    0.7999    0.0012    0.0012    0.0015    0.0015];
GetShapet(xq)

function [] = GetShapet(x)
x_range = linspace(0,1,40);

for i=1:length(x_range)
    h = GetThickness(x_range(i));
    ymin(i) = -h/2;
    yplus(i) = h/2;
end
hold off
plot(x_range, ymin, 'b')
hold on
plot(x_range, yplus, 'b')

%front spar
xc1 = x(1);
ymin1 = -GetThickness(xc1)/2;
yplus1 = GetThickness(xc1)/2;
plot([xc1, xc1], [ymin1, yplus1], 'r')

%aft spar
xc2 = x(2);
ymin2 = -GetThickness(xc2)/2;
yplus2 = GetThickness(xc2)/2;
plot([xc2, xc2], [ymin2, yplus2], 'r')

%aft spar
xcmid = (x(1) + x(2))/2;
yminmid = -GetThickness(xcmid)/2;
yplusmid = GetThickness(xcmid)/2;

% %connect
% xcon = [xc1 xcmid xc2 xc2 xcmid xc1];
% ycon = [yplus1 yplusmid yplus2 ymin2 yminmid ymin1];
% plot(xcon, ycon, 'r')

end



% Initialization
clf, hold off, clear

opt_params;
xc1 = linspace(lb(1), ub(1), 20);
xc2 = linspace(lb(2), ub(2), 20);
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];

%%
for j=1:length(xc1)
  for i=1:length(xc2)
    xq(1) = xc1(j);
    xq(2) = xc2(i);
    [f, b1, b2, g1, g2] = objective(xq);
%     [j, i, xc1(j), xc2(i), f] % for checking
    f_obj(i,j) = f;
    f_g1(i,j) = g1;
    f_g2(i,j) = g2;
    f_b1(i,j) = b1;
    f_b2(i,j) = b2;
    f_tot(i,j) = f + b1 + b2;
  end
end

%%
hold off
% plotting
contour(xc1, xc2, f_obj, 'ShowText', 'on')
xlabel('xc1 [m]'), ylabel('xc2 [m]'), title('objective function contour plot')
grid
hold on

% % g1
contour(xc1, xc2, f_b1, [0.0 0.0], 'b', 'ShowText', 'on')
contour(xc1, xc2, f_b2,  'r', 'ShowText', 'on')
% plot(0.1, 0.6,'r*') %ref design
% % ref lines
% contour(xc1, xc2, f_g1, [0.1 0.1],'b--')
% contour(xc1, xc2, f_g2, [0.01 0.01], 'r--')
% legend('objective function','flutter constraint', 'bending constraint')

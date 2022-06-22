clf, hold off, clear
format long

opt_params;
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];
hi=1e-08;
sq = FiniteDifference(xq, hi);

xc1 = linspace(lb(1), ub(1), 20);
xc2 = linspace(lb(2), ub(2), 20);

for j=1:length(xc1)
  for i=1:length(xc2)
%   Analysis of valve spring.
    xq(1) = xc1(j);
    xq(2) = xc2(i);
    [f, b1, b2, b3] = objective(xq);
    obj_val(j , i) = real(f); 
    barrier1(j , i) = real(b1);
    barrier2(j , i) = real(b2);
    barrier3(j , i) = real(b3);
  end
end

contour(xc1, xc2, obj_val, 'ShowText', 'on')
xlabel('xc1'), ylabel('xc2'), ...
   title('Figure 1     Objective value for k = 0.5')
hold on

%barrier 1
contour(xc1, xc2, barrier1, [0.00 0.00], 'b', 'ShowText', 'on')
contour(xc1, xc2, barrier1, [0.01 0.01], 'b--', 'ShowText', 'on')

%barrier 2
% contour(xc1, xc2, barrier2,  'r', 'ShowText', 'on')
% contour(xc1, xc2, barrier2, [0.01 0.01], 'r--', 'ShowText', 'on')
% 
% %barrier 3
% contour(xc1, xc2, barrier3,  'g', 'ShowText', 'on')
% contour(xc1, xc2, barrier3, [0.05 0.05], 'g--', 'ShowText', 'on')



grid
% Initialization
clear
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

%initial design
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];
% xq = [0.15 0.7 0.0012 0.0012 0.0015 0.0015];

%design starting point
[f, b1, b2, g1, g2]= objective(xq);
obj_old = f;

%beginning function value for ending criteria
tolF = 1e-10;
diff_fval = 1;
fval = 1;
cycle = 0;

% Loop over optimization cycle:
while(diff_fval>tolF && cycle<51)
   	fval_old = fval;
    cycle = cycle + 1;

    % Forward finite diffence gradients of objective function and constraints
    hi=1e-8;        %step size when line searching
    
    % calculating finite difference gradient for the objective function +
    % barriers
    sq = FiniteDifference(xq, hi);
    
    %getting list for plotting
    %
    xc1_lst(cycle) = xq(1);
    xc2_lst(cycle) = xq(2);
    wt1_lst(cycle) = xq(3);
    wt2_lst(cycle) = xq(4);
    ft1_lst(cycle) = xq(5);
    ft2_lst(cycle) = xq(6);
    % Setting of options:
    options = optimset('tolx',1.0e-8,'MaxFunEvals',50);
    
    % Determining bounds. Whichever bound is hit first for a certain alpha
    % is the decider for the limits of the line search
     
    for i=1:length(xq)
        if sq(i) <0
            dif_upper(i) = (xq(i) - lb(i))/abs(sq(i));
            dif_lower(i) = (xq(i) - ub(i))/abs(sq(i));
        else
            dif_upper(i) = (ub(i) - xq(i))/abs(sq(i));
            dif_lower(i) = (lb(i) - xq(i))/abs(sq(i));
        end
    end
    
    %Line search (note the lower and upper bound of alfhaq):
    [alphaq,fval,exitflag] = ...
           fminbnd(@(alpha) objective_mask(alpha, xq,sq), max(dif_lower), 0.5 * min(dif_upper), [options]);
    
    alphaq = alphaq*0.5;
    % Computation of result of line search (new design point):
    for i=1:length(xq)
        xq(i) = xq(i) + alphaq*sq(i);
    end
    
    %function eval
    diff_fval = abs(1 - fval/fval_old);
    [cycle, diff_fval]
    xq_new = xq;
end 
disp('Done')
disp('check constraints, higher than 1 is better')
check_flut = constraints(GetInertia(xq_new))
check_bend = GetInertia(xq_new)

disp('relative value constraints, higher than 1 is better')
relative_flut = check_flut/p_ref
relative_bend = check_bend/(0.8 * I_xx_ref)

disp('design vector')
xq_new

% calculating improvement of design
[f, b1, b2, g1, g2]= objective(xq_new);
obj_new = f;
improvement_percentage  =  (1 - obj_new/obj_old) * 100


%% generating data for spar ranges
xc1 = linspace(lb(1), ub(1), 20);
xc2 = linspace(lb(2), ub(2), 20);
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];

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

hold off
% plotting contour of objective function and constraints
contour(xc1, xc2, f_obj, 'ShowText', 'on') 
hold on
contour(xc1, xc2, f_g1, [0.0 0.0], 'b', 'ShowText', 'on')
contour(xc1, xc2, f_g2, [0.0 0.0], 'k', 'ShowText', 'on')

%plotting path
plot(xc1_lst, xc2_lst, 'r')
plot(xc1_lst(1), xc2_lst(1), 'r', 'Marker', 'o')
plot(xc1_lst(end), xc2_lst(end), 'r', 'Marker', 'o')

%plotting extra contour
contour(xc1, xc2, f_g1, [0.1 0.1],'b--')
contour(xc1, xc2, f_g2, [0.01 0.01],'k--')

%plot settings
xlabel('xc1 [m]'), ylabel('xc2 [m]'), title('objective function contour plot')
legend('objective function', 'flutter constraint', 'bending constraint', 'optimization path')

disp('done')

%% generating web and flange ranges front spar
xc1 = linspace(lb(3), ub(3), 20);
xc2 = linspace(lb(5), ub(5), 20);
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];

for j=1:length(xc1)
  for i=1:length(xc2)
    xq(3) = xc1(j);
    xq(5) = xc2(i);
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

hold off
% plotting contour of objective function and constraints
contour(xc1, xc2, f_obj, 'ShowText', 'on') 
hold on
contour(xc1, xc2, f_g1, [0.0 0.0], 'b', 'ShowText', 'on')
contour(xc1, xc2, f_g2, [0.0 0.0], 'k', 'ShowText', 'on')

%plotting path
plot(wt1_lst, wt2_lst, 'r')
plot(wt1_lst(1), wt2_lst(1), 'r', 'Marker', 'o')
plot(wt1_lst(end), wt2_lst(end), 'r', 'Marker', 'o')

%plotting extra contour
contour(xc1, xc2, f_g1, [0.1 0.1],'b--')
contour(xc1, xc2, f_g2, [0.01 0.01],'k--')

%plot settings
xlabel('ft1 [m]'), ylabel('wt1 [m]'), title('objective function contour plot')
legend('objective function', 'flutter constraint', 'bending constraint', 'optimization path')

disp('done')

%% generating front and aft flange thickness
xc1 = linspace(lb(5), ub(5), 20);
xc2 = linspace(lb(6), ub(6), 20);
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];

for j=1:length(xc1)
  for i=1:length(xc2)
    xq(5) = xc1(j);
    xq(6) = xc2(i);
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

hold off
% plotting contour of objective function and constraints
contour(xc1, xc2, f_obj, 'ShowText', 'on') 
hold on
contour(xc1, xc2, f_g1, [0.0 0.0], 'b', 'ShowText', 'on')
contour(xc1, xc2, f_g2, [0.0 0.0], 'k', 'ShowText', 'on')

%plotting path
plot(ft1_lst, ft2_lst, 'r')
plot(ft1_lst(1), ft2_lst(1), 'r', 'Marker', 'o')
plot(ft1_lst(end), ft2_lst(end), 'r', 'Marker', 'o')

%plotting extra contour
contour(xc1, xc2, f_g1, [0.1 0.1],'b--')
contour(xc1, xc2, f_g2, [0.01 0.01],'k--')

%plot settings
xlabel('ft1 [m]'), ylabel('ft2 [m]'), title('objective function contour plot')
legend('objective function', 'flutter constraint', 'bending constraint', 'optimization path')

disp('done')

%% generating front spar location and flange thickness
xc1 = linspace(lb(1), ub(1), 20);
xc2 = linspace(lb(3), ub(3), 20);
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];

for j=1:length(xc1)
  for i=1:length(xc2)
    xq(1) = xc1(j);
    xq(3) = xc2(i);
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

hold off
% plotting contour of objective function and constraints
contour(xc1, xc2, f_obj, 'ShowText', 'on') 
hold on
contour(xc1, xc2, f_g1, [0.0 0.0], 'b', 'ShowText', 'on')
contour(xc1, xc2, f_g2, [0.0 0.0], 'k', 'ShowText', 'on')

%plotting path
plot(xc1_lst, ft1_lst, 'r')
plot(xc1_lst(1), ft1_lst(1), 'r', 'Marker', 'o')
plot(xc1_lst(end), ft1_lst(end), 'r', 'Marker', 'o')

%plotting extra contour
contour(xc1, xc2, f_g1, [0.1 0.1],'b--')
contour(xc1, xc2, f_g2, [0.01 0.01],'k--')

%plot settings
xlabel('xc1 [m]'), ylabel('ft1 [m]'), title('objective function contour plot')
legend('objective function', 'flutter constraint', 'bending constraint', 'optimization path')

disp('done')




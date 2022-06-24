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

% %initial design
% xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];
% 
% %generating new designs
% lb = [0.05 0.5 0.00075, 0.00075, 0.001, 0.001];
% ub = [0.40 0.8 0.0015 0.0015 0.002 0.002];

des_k = 0.9;
for i=1:length(lb)
    xq(i) = des_k * (ub(i) - lb(i)) + lb(i);
end

%design starting point
obj_old = objective(xq);

%beginning function value for ending criteria
tolF = 1e-10;diff_fval = 1;
fval = 1;
cycle = 0;
% Loop over optimization cycle:
while diff_fval > tolF
    
   	fval_old = fval;
    cycle = cycle + 1
    % Forward finite diffence gradients of objective function and constraints
    hi=1e-8;        %step size when line searching
    
    % calculating finite difference gradient for the objective function +
    % barriers
    sq = FiniteDifference(xq, hi);

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
           fminbnd(@(alpha) objective_mask(alpha, xq,sq), max(dif_lower), min(dif_upper), [options]);
   alphaq
   fval
   exitflag
    
    % Computation of result of line search (new design point):
    for i=1:length(xq)
        xq(i) = xq(i) + alphaq*sq(i);
    end
    
    %function eval
    diff_fval = abs(1 - fval/fval_old);
    xq_new = xq;
end 
disp('Done')
disp('check constraints')
check_flut = constraints(xq_new)

disp('relative value constraints')
relative_flut = check_flut/p_ref

[f, b1, b2] = objective(xq)
obj_new = f + b1 + b2;
improvement_percentage  =  (1 - obj_new/obj_old) * 100


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
objective(xq);

lb = [0.05 0.5 0.00075, 0.00075, 0.001, 0.001];
ub = [0.40 0.8 0.0015 0.0015 0.002 0.002];
obj_old = objective(xq);

%getting bounds for optimizer
lb = zeros(6);
ub = ones(6);

tolF = 0.00001;

%beginning function value for ending criteria
diff_fval = 1;
fval = 1;
cycle = 0;
% Loop over optimization cycle:
while cycle < 10
    
   	fval_old = fval;
    
    cycle = cycle + 1;
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
    
%     dif_min
%     dif_up
    %Line search (note the lower and upper bound of alfhaq):
    [alphaq,fval,exitflag] = ...
           fminbnd(@(alpha) objective_mask(alpha, xq,sq), max(dif_lower), min(dif_upper), [options]);
    
%     % Optimization results:
%     alphaq         % step size
%     fval           % value objective function
    
    
    % Computation of result of line search (new design point):
    for i=1:length(xq)
        xq(i) = xq(i) + alphaq*sq(i);
    end
    
    
    diff_fval = abs(fval - fval_old);
    
end 

check_constraint = constraint2(xq)
relative = check_constraint/p_ref

improvement = objective(xq)/obj_old


opt_params
xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];
hi=1e-08;
sq = FiniteDifference(xq, hi);
for i=1:length(xq)
    if sq(i) <0
        dif_upper(i) = (xq(i) - lb(i))/abs(sq(i));
        dif_lower(i) = (xq(i) - ub(i))/abs(sq(i));
    else
        dif_upper(i) = (ub(i) - xq(i))/abs(sq(i));
        dif_lower(i) = (lb(i) - xq(i))/abs(sq(i));
    end
end

alpha = linspace(max(dif_lower), min(dif_upper), 20);
for i=1:length(alpha)
        % line searching
    x_des = xq + alpha(i)*sq; 

    %objective value
    [f, b1, b2] = objective(x_des);
    f_lst(i) = f;
    g1_lst(i) = g1;
    g2_lst(i) = g2;

end

plot(alpha, g1_lst, 'DisplayName', 'g1')
hold on
plot(alpha, g2_lst, 'DisplayName', 'g2')
legend
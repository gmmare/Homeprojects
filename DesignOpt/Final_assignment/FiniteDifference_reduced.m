function [sq] = FiniteDifference_reduced(x, hi)
% compute finite difference for a given location x and step size hi

%setting change in x
df = zeros(1, 6);

%local objective value
[f, b1, b2] = objective(x);
f_local = f+ b1+ b2;

for i=1:6
    if i<=2
    %setting up new x vector
    dx = zeros(1, 6);
    dx(i) =  hi;

    x_new = x + dx;
    %calculating new objective value
    [f, b1, b2, g1, g2] = objective(x_new);
    f_new = f+ b1+ b2;
    df(i) = (f_new - f_local)/hi;
    end
    
end

sq = -df;
sq = sq/norm(sq) ;

end


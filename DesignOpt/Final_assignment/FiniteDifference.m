function [sq] = FiniteDifference(x, hi)
% compute finite difference for a given location x and step size hi

%setting change in x
df = zeros(1, 6);

%local objective value
f_local = objective(x);

for i=1:length(x)

    %setting up new x vector
    dx = zeros(1, 6);
    dx(i) =  hi;

    x_new = x + dx;
    %calculating new objective value
    f_new = objective(x_new);
    df(i) = (f_new - f_local)/hi;
end

sq = -df;
sq = sq/norm(sq);

end


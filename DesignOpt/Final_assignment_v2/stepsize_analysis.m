%noise analysis

xq = [0.1 0.62 0.0012 0.0012 0.001 0.0015];

%create log space
a = -4;
b = -12;
n = 20;
hi_range = logspace(a,b,n);

df_list = [];
for i= 1:length(hi_range)
    dx = FiniteDifference(xq, hi_range(i));
    df_list = [df_list norm(dx)];
end


semilogx(hi_range, df_list/max(df_list))
ylim([0 1.2])
xlabel('hi [-]'), ylabel('df [-]'), title('Influence of stepsize on gradient')
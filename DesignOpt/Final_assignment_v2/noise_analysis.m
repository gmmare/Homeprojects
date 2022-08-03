%noise analysis

%create 
ixx_range = linspace(1e-8, 1e-6, 20);

flut_list = [];
for i= 1:length(ixx_range)
    flut = constraints(ixx_range(i));
    flut_list = [flut_list flut];
end


semilogx(ixx_range, flut_list/abs(min(flut_list)))
xlabel('ixx'), ylabel('flutter value [-]'), title('Flutter constraint noise analysis')
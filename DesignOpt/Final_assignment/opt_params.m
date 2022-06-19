%parameters used during the optimization
%lower limits
f_spar_l = 0.05;
a_spar_l = 0.51;
f_web_l = 0.0005;
a_web_l = 0.0005;
f_thick_l = 0.0005;
a_thick_l = 0.0005;
bounds_lower = [f_spar_l, a_spar_l, f_web_l, a_web_l, f_thick_l, f_thick_l];

%upper limits
f_spar_u = 0.40;
a_spar_u = 0.80;
f_web_u = 0.0025;
a_web_u = 0.0025;
f_thick_u = 0.0025;
a_thick_u = 0.0025;
bounds_upper = [f_spar_u, a_spar_u, f_web_u, a_web_u, f_thick_u, f_thick_u];

%reference design stiffness and tank crosssectional area
I_xx_ref = 2.8521e-07; %m^4
A_ref = 0.0352; %m^2
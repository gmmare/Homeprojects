%parameters used during the optimization
%{ 
design vector x
1 Front spar location (continuous)
2 Aft spar location (continuous)
3 Web-thickness front spar (continuous)
4 Web-thickness aft spar (continuous)
5 Flange thickness front spar (continuous)
6 Flange thickness aft spar (continuous)
%}

% xq = [0.1 0.6 0.0012 0.0012 0.0015 0.0015];

% flange width:
flange_width = 0.042;

%lower limits
f_spar_l = 0.05;
a_spar_l = 0.35;
f_web_l = 0.0005;
a_web_l = 0.0005;
f_thick_l = 0.0005;
a_thick_l = 0.0005;
lb = [f_spar_l, a_spar_l, f_web_l, a_web_l, f_thick_l, f_thick_l];

%upper limits
f_spar_u = 0.25;
a_spar_u = 0.80;
f_web_u = 0.0025;
a_web_u = 0.0025;
f_thick_u = 0.0025;
a_thick_u = 0.0025;
ub = [f_spar_u, a_spar_u, f_web_u, a_web_u, f_thick_u, f_thick_u];

%reference design stiffness and tank crosssectional area
I_xx_ref = 2.8521e-07; %m^4
A_ref = 0.035; %m^2
p_ref = -0.003528; %[-]
W_ref = 3.9968e-04; % [m^2]based on crossection of structural parts

% objective weight:
k = 0.5;
c = 0.25; %correction factor weights
r = 800;


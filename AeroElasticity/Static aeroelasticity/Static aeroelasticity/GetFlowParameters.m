function flow = GetFlowParameters

rho       = 1.225;
V         = 23;
alpha0    = 5*pi/180;
C_M_AC    = 0;
C_L_alpha = 2*pi;

q = .5*rho*V^2;

flow.rho       = rho;
flow.V         = V;
flow.alpha0    = alpha0;
flow.C_M_AC    = C_M_AC;
flow.C_L_alpha = C_L_alpha;
flow.q         = q;
function [Ka,F_st_alpha0,F_st_M_AC] = GetSteadyFlowForces(TS,FL)

% Aerodynamic stiffness per unit q
Ka = [0      -TS.S*FL.C_L_alpha                  ;
      0      TS.S*FL.C_L_alpha*(.5+TS.a)*TS.b    ];

% Aero forces per unit q per unit alpha_0
F_st_alpha0 = [-TS.S*FL.C_L_alpha;TS.S*FL.C_L_alpha*(.5+TS.a)*TS.b];

% Aero moment per unit q
F_st_M_AC   = [0;TS.S*FL.C_M_AC*2*TS.b];
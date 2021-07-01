function [Ka,Ca,F_st_alpha0,F_st_M_AC] = GetQuasiSteadyFlowForces(TS,FL)

% Aerodynamic stiffness per unit q
Ka = [0      -TS.S*FL.C_L_alpha                  ;
      0      TS.S*FL.C_L_alpha*(.5+TS.a)*TS.b    ];
  
% Aerodynamic damping per unit q/V
Ca = [-TS.S*FL.C_L_alpha                     0    ;
      TS.S*FL.C_L_alpha*(0.5*TS.a)*TS.b      0    ];

% Aero forces per unit q per unit alpha_0
F_st_alpha0 = [-TS.S*FL.C_L_alpha;TS.S*FL.C_L_alpha*(.5+TS.a)*TS.b];

% Aero moment per unit q
F_st_M_AC   = [0;TS.S*FL.C_M_AC*2*TS.b];
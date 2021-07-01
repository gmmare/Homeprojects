function varargout = GetTypicalSectionStructuralMatrices(parameters)

Ms = [   parameters.m                     parameters.S_theta                                                       parameters.S_beta                                       ;
      parameters.S_theta                  parameters.I_theta                                       parameters.I_beta+(parameters.c-parameters.a)*parameters.b*parameters.S_beta    ;
      parameters.S_beta    parameters.I_beta+(parameters.c-parameters.a)*parameters.b*parameters.S_beta                    parameters.I_beta                                       ];

Ks = [parameters.K_h         0                    0          ;
        0           parameters.K_theta            0          ;
        0                  0            parameters.K_beta    ];
    
varargout{1} = Ms;
varargout{2} = Ks;
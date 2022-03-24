clear all
close all
clc

flow  = GetFlowParameters;
param = GetTypicalSectionParameters;
Ms = [param.m param.S_theta;param.S_theta param.I_theta];
Ks = [param.K_h 0;0 param.K_theta];

kvec = linspace(.1,5,100);

out = [];
for i=1:length(kvec)
    k = kvec(i);
    A = GetTheodorsenMatrix(param,k);
    B = Ks\(Ms+1/2*flow.rho*param.b^2*A/k^2);
    [~,D] = eig(B);
    lambda = diag(D);
    ome = (1./real(lambda)).^0.5;
    g   = ome.^2.*imag(lambda);
    V   = ome*param.b/k;
    out = [out;k  ome(1)  g(1)  V(1)  ome(2)  g(2)  V(2)];
end

plot(out(:,4),out(:,3)), hold on
plot(out(:,7),out(:,6))
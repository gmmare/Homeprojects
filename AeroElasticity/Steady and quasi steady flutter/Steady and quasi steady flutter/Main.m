clear all
close all
clc

%% Get flow parameters
FL = GetFlowParameters;

%% Get typical section parameters
TS = GetTypicalSectionParameters;

%% Build structural matrices
[Ms,Ks] = GetTypicalSectionStructuralMatrices(TS);

%% Build steady aerodynamic model
% [Ka,F_st_alpha0,F_st_M_AC] = GetSteadyFlowForces(TS,FL);

% %% Build quasi steady aerodynamic model
[Ka,Ca,F_st_alpha0,F_st_M_AC] = GetQuasiSteadyFlowForces(TS,FL);

%% Calculate eigenfrequencies with changing flow velocity
q = linspace(eps,150,100);
V = sqrt(2*q/FL.rho);
syms ps
p_save = [];
for i=1:length(q)
    Kae = Ks(1:2,1:2)-q(i)*Ka;
%     A = ps^2*Ms(1:2,1:2)+Kae;
    A = ps^2*Ms(1:2,1:2)-ps*q(i)/V(i)*Ca+Kae;
    DA = det(A);
    characteristic_equation = sym2poly(DA);
    p = roots(characteristic_equation);
    p_save = [p_save p];
end

%% Plot flutter diagrams
figure, hold on
plot(q,imag(p_save),'.','LineWidth',2)
xlabel('Dynamic pressure [N/m^2]')
ylabel('Frequency \omega [1/s]')
yl = ylim;
axis([xlim 0 yl(2)])
set(gca,'FontSize',14)

figure, hold on
plot(q,real(p_save),'.','LineWidth',2)
xlabel('Dynamic pressure [N/m^2]')
ylabel('Damping \sigma [1/s]')
plot(xlim,[0 0],'k','LineWidth',2)
set(gca,'FontSize',14)

figure, hold on
plot(real(p_save),imag(p_save),'r.','LineWidth',2)
xlabel('Damping \sigma [1/s]')
ylabel('Frequency \omega [1/s]')
plot([0 0],ylim,'k','LineWidth',2)
set(gca,'FontSize',14)
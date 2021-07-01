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
[Ka,F_st_alpha0,F_st_M_AC] = GetSteadyFlowForces(TS,FL);

%% Iterative aeroelastic solution
tol    = 1;
iter   = 0;
x_old  = [eps;eps];
x_save = [];
while(tol>1e-6)
    iter = iter+1;
    x_new = Ks(1:2,1:2)\(FL.q*F_st_alpha0*(FL.alpha0+x_old(2))+FL.q*F_st_M_AC);
    tol   = abs(1-norm(x_new)/norm(x_old));
    x_old = x_new;
    x_save = [x_save x_old];
end

%% Monolithic aeroelastic solution
x_mono = (Ks(1:2,1:2)-FL.q*Ka)\(FL.q*F_st_alpha0*FL.alpha0+FL.q*F_st_M_AC);

%% Plot convergence history
figure, hold on
plot(1:iter,x_save(1,:),'LineWidth',2)
plot([0 iter],[x_mono(1,1) x_mono(1,1)],'LineWidth',2)
xlabel('Number of iterations','FontSize',14)
ylabel('Heave displacement [m]','FontSize',14)
legend('Partitioned','Monolithic')
set(gca,'Fontsize',14)
axis([0 iter -5 0])

figure, hold on
plot(1:iter,x_save(2,:)*180/pi,'LineWidth',2)
plot([0 iter],180/pi*[x_mono(2,1) x_mono(2,1)],'LineWidth',2)
legend('Partitioned','Monolithic','Location','SouthEast')
xlabel('Number of iterations','FontSize',14)
ylabel('Twist rotation [deg]','FontSize',14)
set(gca,'Fontsize',14)
axis([0 iter 0 7])
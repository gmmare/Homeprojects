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

%% Iterative aeroelastic trim solution
tol1        = 1;
iter1       = 0;
alpha0_new  = 1*pi/180;
alpha0_old  = alpha0_new;
x_old       = [eps;eps];
iter2_save  = [];
alpha0_save = alpha0_old;
x_save      = x_old;
while(tol1>1e-6)    
    tol2    = 1;
    iter2   = 0;
    iter1   = iter1+1;
    while(tol2>1e-6)
        iter2 = iter2+1;
        x_new = Ks(1:2,1:2)\(FL.q*F_st_alpha0*(alpha0_new+x_old(2))+FL.q*F_st_M_AC);
        tol2  = abs(1-norm(x_new)/norm(x_old));
        x_old = x_new;
        x_save = [x_save x_old];
    end
    iter2_save = [iter2_save iter2];
    Lift = FL.q*TS.S*FL.C_L_alpha*(x_new(2)+alpha0_old);
    alpha0_new = TS.W/Lift*alpha0_old;
    tol1 = abs(1-alpha0_new/alpha0_old);
    alpha0_save = [alpha0_save alpha0_new];
    alpha0_old = alpha0_new;
end

%% Monolithic aeroelastic solution
x_mono = [(Ks(1:2,1:2)-FL.q*Ka) [FL.q*TS.S*FL.C_L_alpha;-FL.q*TS.S*FL.C_L_alpha*(.5+TS.a)*TS.b];[0 FL.q*TS.S*FL.C_L_alpha FL.q*TS.S*FL.C_L_alpha]]\[0;FL.q*F_st_M_AC(2);TS.W];

%% Plot convergence history
figure, hold on
plot(x_save(1,:),'LineWidth',2)
plot([0 length(x_save(1,:))],[x_mono(1,1) x_mono(1,1)],'LineWidth',2)
legend('Partitioned','Monolithic','Location','SouthEast')
xlabel('Number of iterations','FontSize',14)
ylabel('Heave displacement [m]','FontSize',14)
set(gca,'Fontsize',14)
axis([0 length(x_save(1,:)) -5 0])

figure, hold on
plot(x_save(2,:)*180/pi,'LineWidth',2)
plot([0 length(x_save(2,:))],180/pi*[x_mono(2,1) x_mono(2,1)],'LineWidth',2)
legend('Partitioned','Monolithic','Location','SouthEast')
xlabel('Number of iterations','FontSize',14)
ylabel('Twist rotation [deg]','FontSize',14)
set(gca,'Fontsize',14)
axis([0 length(x_save(2,:)) 0 7])

figure, hold on
plot(alpha0_save*180/pi,'LineWidth',2)
plot([0 length(alpha0_save)],180/pi*[x_mono(3,1) x_mono(3,1)],'LineWidth',2)
legend('Partitioned','Monolithic','Location','SouthEast')
xlabel('Number of iterations','FontSize',14)
ylabel('Rigid angle of attack [deg]','FontSize',14)
set(gca,'Fontsize',14)
axis([0 length(alpha0_save) 0 (max(alpha0_save))*180/pi+1])
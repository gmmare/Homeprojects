clear all
close all
clc

%% Finite element settings

elem = [];
conv = [];
div_list = [];
for k=1:40
elem = [elem, k];

Ne = k;            % Number of elements
Nn = Ne+1;          % Number of nodes
Ndof = 3*(Ne+1);    % Number of degrees of freedom
% The fixed degrees of freedom assume that you arranged the DOFs as:
% [displacements | velocities | aerodynamic states]. If you have another
% arrangement in your mode, please update the following line to reflect
% your own situation.
fxdof = [1:3 3*Nn+(1:3) 6*Nn+(1:4)]; % Constrained degrees of freedom
frdof = setdiff(1:Ndof,fxdof);   % Unconstrained degrees of freedom

Ks  = zeros(Ndof,Ndof);
Ms  = zeros(Ndof,Ndof);
EFT = zeros(Ne,6);

for i=1:Ne
    EFT(i,1:6) = [1:6]+(i-1)*3; % Element freedom table, element connectivity
end

%% Stick model input
span   = 16; %half span
L      = span/Ne;
E      = 70e9;
G      = E/2/(1+0.3);
Ixx    = 2e4/E;
J      = 1e4/G;

Itheta = 0.1;
m      = 0.75;

c          = 1;
S          = zeros(Nn,1);
S(2:end-1) = L*c;
S(1)       = L*c/2;
S(end)     = L*c/2;
a          = 0;
b          = c/2;
xtheta     = 0;

rho = 0.0889;

%% Structural stiffness matrix
Ke = zeros(6,6);
Kb = [  (12*E*Ixx)/L^3,  (6*E*Ixx)/L^2, -(12*E*Ixx)/L^3,  (6*E*Ixx)/L^2;
    (6*E*Ixx)/L^2,    (4*E*Ixx)/L,  -(6*E*Ixx)/L^2,    (2*E*Ixx)/L;
    -(12*E*Ixx)/L^3, -(6*E*Ixx)/L^2,  (12*E*Ixx)/L^3, -(6*E*Ixx)/L^2;
    (6*E*Ixx)/L^2,    (2*E*Ixx)/L,  -(6*E*Ixx)/L^2,    (4*E*Ixx)/L];
Kt = [  (G*J)/L, -(G*J)/L;
    -(G*J)/L,  (G*J)/L];
Ke([1:2,4:5],[1:2,4:5]) = Kb;
Ke([3,6],[3,6])         = Kt;

for i=1:Ne    
    Ks(EFT(i,:),EFT(i,:)) = Ks(EFT(i,:),EFT(i,:))+Ke;   % Assembly of the global stiffness matrix
end

%% Structural mass matrix
Me = zeros(6,6);
Mbb = [     (13*L*m)/35, (11*L^2*m)/210,      (9*L*m)/70, -(13*L^2*m)/420;
    (11*L^2*m)/210,    (L^3*m)/105,  (13*L^2*m)/420,    -(L^3*m)/140;
    (9*L*m)/70, (13*L^2*m)/420,     (13*L*m)/35, -(11*L^2*m)/210;
    -(13*L^2*m)/420,   -(L^3*m)/140, -(11*L^2*m)/210,     (L^3*m)/105];
Mtb = [ (7*L*m*xtheta)/20, (L^2*m*xtheta)/20, (3*L*m*xtheta)/20, -(L^2*m*xtheta)/30;
    (3*L*m*xtheta)/20, (L^2*m*xtheta)/30, (7*L*m*xtheta)/20, -(L^2*m*xtheta)/20];
Mbt = [  (7*L*m*xtheta)/20,  (3*L*m*xtheta)/20;
    (L^2*m*xtheta)/20,  (L^2*m*xtheta)/30;
    (3*L*m*xtheta)/20,  (7*L*m*xtheta)/20;
    -(L^2*m*xtheta)/30, -(L^2*m*xtheta)/20];
Mtt = [ (Itheta*L)/3, (Itheta*L)/6;
    (Itheta*L)/6, (Itheta*L)/3];
Me([1:2,4:5],[1:2,4:5]) = Mbb;
Me([3,6],[3,6])         = Mtt;
Me([3,6],[1:2,4:5])     = Mtb;
Me([1:2,4:5],[3,6])     = Mbt;

for i=1:Ne
    Ms(EFT(i,:),EFT(i,:)) = Ms(EFT(i,:),EFT(i,:))+Me;
end

%% structural eigenfrequencies
[V, D] = eig(Ks(4:end, 4:end), Ms(4:end, 4:end));
coupled = sqrt(diag(D)); %eigenfrequencies in hz, order is heave, bending twist
conv = [conv, coupled(1:3)];



%% divergence speed
n_dof = (Ne+1) * 3;
Ka_steady = zeros(Ne+1,Ne+1); %Steady state aerodynamic matrix
for i=1:Ne + 1
    Ka_steady(1 + 3*(i - 1), 3 + 3*(i-1)) = 2*pi * c * L; %lift 
    Ka_steady(3 + 3*(i - 1), 3 + 3*(i-1)) = 2*pi * c *(0.5 +a)*b * L; %moment
end
Ka_steady(end -2, end) = 0.5 * 2*pi * c * L; %lift 
Ka_steady(end, end) = 0.5 * 2*pi * c *(0.5 +a)*b * L; %moment
qlist = eig(Ks(4:end, 4:end), Ka_steady(4:end, 4:end));
qlist = qlist(qlist>0.01 & isfinite(qlist));

v_div = sqrt((2*min(qlist))/rho);
div_list = [div_list, v_div];
end

figure(1)
%plotting struc eigenfunction convergence
plot(elem, conv)
legend('1st ','2nd', '3rd')
ylabel('Frequency[rad/s]') 
xlabel('Number of wing elements [-]') 
hold off

figure(2)
%plotting divergence speed convergence
text = join('Converged result: ' + string(div_list(end)) + ' [ms^{-1}]');
plot(elem, div_list)
yline(div_list(end),'-', text)
ylabel('Divergence speed [ms^-1]') 
xlabel('Number of wing elements [-]') 

%% plotting divergence mode
lamda = 0.5 * rho * v_div^2;
M = Ks(4:end, 4:end) - lamda*Ka_steady(4:end, 4:end);
v_zero = zeros(elem(end)*3, 1);
x = null(M);
heave = [];
bending = [];
twist = [];
for i=1:length(elem)
    heave = [heave, -x(1 + 3*(i - 1))];
    bending = [bending, -x(2 + 3*(i - 1))];
    twist = [twist, -x(3 + 3*(i - 1))];
end
% 
% heave = heave/max(heave);
% bending = bending/max(bending);
% twist = twist/max(twist);
elem = elem/max(elem);

figure(3)
plot(elem, heave,'DisplayName', 'heave')
hold on
plot(elem, bending,'DisplayName', 'bending')
plot(elem, twist,'DisplayName', 'twist')

ylabel('relative displacement[-]') 
xlabel('spanwise position [-]') 
legend

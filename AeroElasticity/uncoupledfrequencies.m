clear all
close all
clc

%% Finite element settings

freq_list = [];

xtheta_list = linspace(-0.35, 0.35, 100);
for i=1:length(xtheta_list)
xtheta     = xtheta_list(i);
Ne = 12;            % Number of elements
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
coupled(1:3);
freq_list = [freq_list, coupled(1:3)];
end

figure(1)
%plotting struc eigenfunction convergence
plot(xtheta_list, freq_list)
legend('Heave','bending', 'twist')
ylabel('Frequency[rad/s]') 
xlabel('xtheta [-]')







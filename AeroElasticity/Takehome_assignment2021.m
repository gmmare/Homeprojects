clear all
close all
clc

%% Finite element settings
Ne = 16;            % Number of elements
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

%% Steady aerodynamics, calculating trimming
Vinf = 25;
q = 0.5 * rho * Vinf^2 ; % Dynamic pressure, vary this value to assess aeroelastic instability
V_range = [];
alpha_range = [];
alpha_range_rigid = [];
for k=20:0.1:25
q = 0.5 * rho * k^2;
V_range = [V_range,k]; 

%non rigid
A = [];
n_dof = (Ne+1) * 3;
Ka_steady = zeros(Ne+1,Ne+1); %Steady state aerodynamic matrix
for i=1:Ne + 1
    Ka_steady(1 + 3*(i - 1), 3 + 3*(i-1)) = q * 2*pi * c * L; %lift 
    Ka_steady(3 + 3*(i - 1), 3 + 3*(i-1)) = q * 2*pi * c *(0.5 +a)*b * L; %moment
end
Ka_steady(end -2, end) = 0.5 * q * 2*pi * c * L; %lift 
Ka_steady(end, end) = 0.5 * q * 2*pi * c *(0.5 +a)*b * L; %moment

% difference stiffness matrices
K_dif = Ks - Ka_steady;

% extra conditions due to alpha 0 per section
section_terms = zeros(n_dof, 1);
for i=1:Ne + 1
    section_terms(1 + 3*(i - 1)) = -q * 2*pi * c * L;    %lift term
    section_terms(3 + 3*(i - 1)) = -q * 2*pi * c *(0.5 + a)*b * L;    %moment term
end

%rigid
A_r = [];
Ka_rigid = zeros((Ne+1)*3,(Ne+1)*3); %Steady state aerodynamic matrix
% difference stiffness matrices

K_dif_rigid = Ks - Ka_rigid;

% extra conditions due to alpha 0 per section
section_terms = zeros(n_dof, 1);
for i=1:Ne + 1
    section_terms(1 + 3*(i - 1)) = -q * 2*pi * c * L;    %lift term
    section_terms(3 + 3*(i - 1)) = -q * 2*pi * c *(0.5 + a)*b * L;    %moment term
end


% extra equation for trim condition
trim_terms = zeros(1, n_dof +1);
trim_terms_rigid = zeros(1, n_dof +1);
for i=1:Ne + 1
    trim_terms(3 + 3*(i - 1)) = q * 2*pi * c * L;    %lift term, for dof, positive downwards
end
trim_terms(end-1) = q * 2*pi * c * L*0.5;
trim_terms(end) = q * 2*pi*span*c;
trim_terms_rigid(end) = q * 2*pi*span*c;

% difference stiffness matrices
K_diff = Ks - Ka_steady;

% getting the total matrix for trim
K_trim = cat(2, K_diff, section_terms);
K_trim = cat(1, K_trim, trim_terms);
K_trim_rigid = cat(2, K_dif_rigid, section_terms);
K_trim_rigid = cat(1, K_trim_rigid, trim_terms_rigid);

k_trim_final = K_trim(4:end, 4:end);
k_trim_final_rigid = K_trim_rigid(4:end, 4:end);


% weight vector
v_weight = zeros(n_dof-3 + 1, 1);
v_weight(end) = 730/2;

x_solve = linsolve(k_trim_final, v_weight);
x_solve_rigid = linsolve(k_trim_final_rigid, v_weight);
alpha_0 = x_solve(end) * 180 /pi;
alpha_0_rigid = x_solve_rigid(end) * 180 /pi;

alpha_range = [alpha_range, alpha_0];
alpha_range_rigid = [alpha_range_rigid, alpha_0_rigid];
% % plotting wing deflection of trimmed state
% deflection = zeros(Ne);
% section = zeros(Ne);
% for i=1:Ne
%     deflection(i) = x_solve(3 + 3*(i - 1),1); %positive downwards 
%     section(i) = i;
% end

end
plot(V_range, alpha_range, 'r-o'); hold on
plot(V_range, alpha_range_rigid, 'b-+');
legend('Flexibel','Rigid')
xticks([20, 21, 22, 23, 24, 25])
xlabel('Flightspeed freesteam velocity [ms^-1]') 
ylabel('Trim angle of attack [degrees]') 

%% Reduce the system by applying the boundary conditions
% Make sure you order the DOFs such that the first three are the clamped
% displacement, bending rotation and torsional rotation
Ar = A;
% Ar(fxdof,fxdof) = [];

%% Calculate the eigenvalues, evaluate and verify
[Vr,Dr] = eig(Ar);
Dr      = diag(Dr);
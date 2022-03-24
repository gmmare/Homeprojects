clear all
close all
clc

%% Finite element settings
Ne = 30;            % Number of elements
Nn = Ne+1;          % Number of nodes
Ndof = 3*(Ne+1);    % Number of degrees of freedom
% The fixed degrees of freedom assume that you arranged the DOFs as:
% [displacements | velocities | aerodynamic states]. If you have another
% arrangement in your mode, please update the following line to reflect
% your own situation.
fxdof = [1:3 3*Nn+(1:3) 6*Nn+(1:4)]; % Constrained degrees of freedom
frdof = setdiff(1:Ndof,fxdof);   % Unconstrained degrees of freedom

V = 24 ;
rho = 0.0889;
q = 0.5 * rho * V^2;

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

%Wagner 
psi1_w = 0.165;
psi2_w = 0.335;
eps1_w = 0.0455;
eps2_w = 0.3;
psi1_k = 0.5;
psi2_k = 0.5;
eps1_k = 0.13;
eps2_k = 1;

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
% Building aerodynamic forces into state space matrices
n_dof = (Ne+1) * 3;

% aerodynamic mass
Ma = zeros(n_dof,n_dof);
for i=1:Ne + 1
    Ma(1 + 3*(i - 1), 1 + 3*(i-1)) = -1;  
    Ma(1 + 3*(i - 1), 3 + 3*(i-1)) = -b*a; 
    Ma(3 + 3*(i - 1), 1 + 3*(i-1)) = -b*a; 
    Ma(3 + 3*(i - 1), 3 + 3*(i-1)) = -b^2 *(1/8 + a^2); 
end
Ma = Ma(4:end, 4:end) * 0.5 * 2 * pi * rho * b^2 * L; %remove root node
Ma(end-2:end, end-2:end) = Ma(end-2:end, end-2:end) * 0.5;

% aerodynamic damping
Ca1 = zeros(n_dof,n_dof);
for i=1:Ne + 1
    Ca1(1 + 3*(i - 1), 3 + 3*(i-1)) = V; 
    Ca1(3 + 3*(i - 1), 3 + 3*(i-1)) = -V*b*(1/2 + a); 
end

Ca2 = zeros(n_dof,n_dof);
for i=1:Ne + 1
    Ca2(1 + 3*(i - 1), 1 + 3*(i-1)) = -(1-psi1_w-psi2_w)/V;
    Ca2(1 + 3*(i - 1), 3 + 3*(i-1)) = b*(1/2 - a)*(1-psi1_w-psi2_w)/V; 
    Ca2(3 + 3*(i - 1), 1 + 3*(i-1)) = -b*(1/2 + a)*(1-psi1_w-psi2_w)/V; 
    Ca2(3 + 3*(i - 1), 3 + 3*(i-1)) = b*(1/2 + a)*b*(1/2 - a)*(1-psi1_w-psi2_w)/V; 
end
Ca1 = Ca1(4:end, 4:end)* 0.5 * 2 * pi * rho * b^2 * L;%remove root node
Ca2 = Ca2(4:end, 4:end)* 2 * 2 * pi * q * b * L;%remove root node
Ca = Ca1 + Ca2;
Ca(end-2:end, end-2:end) = Ca(end-2:end, end-2:end) * 0.5; %half node at end
% aerodynamic stiffness
Ka = zeros(n_dof,n_dof); %Steady state aerodynamic matrix
for i=1:Ne + 1
    Ka(1 + 3*(i - 1), 1 + 3*(i-1)) = -(psi1_w * eps1_w/b + psi2_w*eps2_w/b); 
    Ka(1 + 3*(i - 1), 3 + 3*(i-1)) = (psi1_w * eps1_w + psi2_w*eps2_w)*(1/2 - a) +(1-psi1_w-psi2_w);
    Ka(3 + 3*(i - 1), 1 + 3*(i-1)) = -(psi1_w * eps1_w/b + psi2_w*eps2_w/b)*b*(1/2 + a);
    Ka(3 + 3*(i - 1), 3 + 3*(i-1)) = ((psi1_w * eps1_w + psi2_w*eps2_w)*(1/2 - a)+(1-psi1_w-psi2_w))*b*(1/2 + a);
end
Ka = Ka(4:end, 4:end)* 2 * 2 * pi * q * b * L;%remove root node
Ka(end-2:end, end-2:end) = Ka(end-2:end, end-2:end) * 0.5;

% lag terms
W = zeros(n_dof,(Ne+1) * 4);
for i=1:Ne + 1
    W(1 + 3*(i - 1), 1 + 4*(i-1)) = psi1_w * eps1_w^2 * V/(b^2);
    W(1 + 3*(i - 1), 2 + 4*(i-1)) = psi2_w * eps2_w^2 * V/(b^2); 
    W(1 + 3*(i - 1), 3 + 4*(i-1)) = psi1_w * eps1_w * V/b * (1 - eps1_w*(0.5 - a));
    W(1 + 3*(i - 1), 4 + 4*(i-1)) = psi2_w * eps2_w * V/b * (1 - eps2_w*(0.5 - a));
    W(3 + 3*(i - 1), 1 + 4*(i-1)) = psi1_w * eps1_w^2 * V/(b^2) * b*(1/2 + a);
    W(3 + 3*(i - 1), 2 + 4*(i-1)) = psi2_w * eps2_w^2 * V/(b^2) * b*(1/2 + a); 
    W(3 + 3*(i - 1), 3 + 4*(i-1)) = psi1_w * eps1_w * V/b * (1 - eps1_w*(0.5 - a)) * b*(1/2 + a);
    W(3 + 3*(i - 1), 4 + 4*(i-1)) = psi2_w * eps2_w * V/b * (1 - eps2_w*(0.5 - a)) * b*(1/2 + a);
end
W = W(4:end, 5:end)* 2 * 2 * pi * q * b * L;%remove root node
W(end-2:end, end-3:end) = W(end-2:end, end-3:end) * 0.5;

%gust lag terms
Wgust = zeros(n_dof,2);
for i=1:Ne + 1
    Wgust(1 + 3*(i - 1), 1 ) = psi1_k * eps1_k;
    Wgust(1 + 3*(i - 1), 2 ) = psi2_k * eps2_k; 
    Wgust(3 + 3*(i - 1), 1 ) = psi1_k * eps1_k * b*(1/2 + a);
    Wgust(3 + 3*(i - 1), 2 ) = psi2_k * eps2_k * b*(1/2 + a); 

end
Wgust = Wgust(4:end,:)* 2 * 2 * pi * q * b * L;%remove root node
Wgust(end-2:end,:) = Wgust(end-2:end,:) * 0.5;

% lag term derivative Q 
Q = zeros((Ne+1) * 4,(Ne+1) * 4);
for i=1:Ne + 1
    Q(1 + 4*(i - 1), 1 + 4*(i-1)) = -eps1_w * V/b;
    Q(2 + 4*(i - 1), 2 + 4*(i-1)) = -eps2_w * V/b; 
    Q(3 + 4*(i - 1), 3 + 4*(i-1)) = -eps1_w * V/b;
    Q(4 + 4*(i - 1), 4 + 4*(i-1)) = -eps2_w * V/b;
end
Q = Q(5:end, 5:end);%remove root node

% lag term derivative R 
R = zeros(2,2);
for i=1:Ne + 1
    Q(1 , 1) = -eps1_k * V/b;
    Q(2 , 2) = -eps2_k * V/b; 
end

% lag term derivative P 
P = zeros((Ne+1) * 4,n_dof);
for i=1:Ne + 1
    P(1 + 4*(i - 1), 1 + 3*(i-1)) = 1;
    P(2 + 4*(i - 1), 1 + 3*(i-1)) = 1; 
    P(3 + 4*(i - 1), 3 + 3*(i-1)) = 1;
    P(4 + 4*(i - 1), 3 + 3*(i-1)) = 1;
end
P = P(5:end, 4:end);%remove root node

% lag term derivative S 
S = ones(2,1);

% extra matrices
I = eye(Ne * 3);
zero_22 = zeros(Ne *3);
zero_23 = zeros(Ne*3, Ne*4);
zero_24 = zeros(Ne*3, 2);
zero_31 = zeros(Ne * 4, Ne *3);
zero_34 = zeros(Ne * 4, 2);
zero_41 = zeros(2, Ne*3);
zero_42 = zeros(2, Ne*3);
zero_43 = zeros(2, Ne*4);


%% building the matrices
Mae = Ms(4:end, 4:end) - Ma;
Kae = Ks(4:end, 4:end) - Ka;

%first row
A1 = cat(2, Mae\Ca, -Mae\Kae);
A1 = cat(2, A1, Mae\W);
A1 = cat(2, A1, Mae\Wgust);
%2nd row
A2 = cat(2, I, zero_22);
A2 = cat(2, A2, zero_23);
A2 = cat(2, A2, zero_24);
%3rd row
A3 = cat(2, zero_31, P);
A3 = cat(2, A3, Q);
A3 = cat(2, A3, zero_34);
%4th row
A4 = cat(2, zero_41, zero_42);
A4 = cat(2, A4, zero_43);
A4 = cat(2, A4, R);

% Assembling total matrix
A = cat(1, A1, A2);
A = cat(1, A, A3);
A = cat(1, A, A4);

%% building B matrix
Fgust = zeros(Ne * 3, 1);
for i=1:Ne
    Fgust(1 + 3*(i-1)) = 1 - psi1_k - psi2_k;
    Fgust(3 + 3*(i-1)) = (1 - psi1_k - psi2_k) * b * (0.5 + a);
end 
Fgust = Fgust * 2 * 2 * pi * q * b * L;
Fgust(end-2:end) = Fgust(end-2:end) * 0.5;
B1 = Mae\Fgust;

zero_21 = zeros(Ne * 3, 1);
zero_31 = zeros(Ne * 4, 1);


B = cat(1, B1, zero_21);
B = cat(1, B, zero_31);
B = cat(1, B, S);

C = zeros(Ne * 10 + 2, Ne * 10 + 2);
D = zeros(Ne * 10 + 2, 1);

%% building state space
sys  = ss(A, B, C, D);


%% time response
Wg = 0.3;
response_time = 10;
gust_length = 3; %seconds
t_gust = linspace(0,gust_length, gust_length*50);
u1 = 1 - cos(t_gust * 50 *2* pi/(length(t_gust)));
u = cat(2, u1, zeros(1, response_time*50));
t = linspace(0,gust_length + response_time, (gust_length + response_time)*50);
u2 = u * Wg/V;
[y, t, x] = lsim(sys, u2, t);

% u = rad2deg(atan(u));
figure(1)
yyaxis right
xlabel('time [t]')
plot(t, u) 
ylabel('Gust input [degr]')

yyaxis left
plot(t, x(:,88))
ylabel('tip displacement [m]')


%% gust length variation
gust_plots = [];
gust_lengths = linspace(1,8,40);
response_time = 40;
gust_amp = [0.2, 0.4, 0.6];
for i=1:length(gust_amp)
    Wg = gust_amp(i);
    length_disp = [];
    for i=1:length(gust_lengths)
    gust_length = gust_lengths(i); %seconds
    t_gust = linspace(0,gust_length, gust_length*50);

    u = 1 - cos(t_gust * 50* 2 * pi/(length(t_gust)));
    u = cat(2, u, zeros(1, response_time*50));
    u = u * Wg/V;
    t = linspace(0,gust_length + response_time, (gust_length + response_time)*50);
    [y, t, x] = lsim(sys, u, t);
    length_disp = [length_disp, max(x(:,88))];
    end
    gust_plots = [gust_plots; length_disp];

end

figure(2)
plot(gust_lengths * V, gust_plots)
xlabel('Gust length [m]')
ylabel('Maximum tip displacement [m]')
yline(1)
legend('Wg=0.2', 'Wg=0.4', 'Wg=0.6')




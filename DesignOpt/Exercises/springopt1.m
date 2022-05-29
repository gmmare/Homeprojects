% Initialization
clf, hold off, clear

% Assignment of constant design parameter values
springparams1;

% max params
tau12max = 600*10^6;
ref_mass = 0.04;

%% Matrix of output values for combinations of design variables D and d 
D = [0.020:0.001:0.040];
d = [0.002:0.0002:0.005];
for j=1:1:length(d)
  for i=1:1:length(D)
%   Analysis of valve spring.
    [svol,smass,bvol,matc,manc,Lmin,L2,k,F1,F2,Tau1,Tau2,freq1]=...
    springanalysis1(D(i),d(j),L0,L1,n,E,G,rho,Dv,h,p1,p2,nm,ncamfac,nne,matp,bldp);
    spring_mass = obj_function(smass)
    [t1, f1, f2, L_spring, freq] = constraint(Tau2, F1, F2, L2, Lmin, freq1, Dv, p1, p2, nm)
    funk_mass(j,i) = spring_mass;
    funk_tau(j,i) = t1;
    funk_f1(j,i) = f1;
    funk_f2(j,i) = f2;
    funk_length(j,i) = L_spring;
    funk_freq(j,i) = freq;
  end
end

% plotting
contour(D, d, funk_mass, [ 0.03 0.04 0.05 0.06], 'ShowText', 'on')
xlabel('D (m)'), ylabel('d (m)'), title('Spring mass')
grid
hold on
% tau
contour(D, d, funk_tau, [0.0 0.0], 'b')
contour(D, d, funk_tau, [0.01 0.01], 'b--')
% forces
contour(D, d, funk_f1, [0.0 0.0], 'r')
contour(D, d, funk_f1, [0.01 0.01], 'r--')
contour(D, d, funk_f2, [0.0 0.0], 'g')
contour(D, d, funk_f2, [0.01 0.01], 'g--')
% length
contour(D, d, funk_length, [0.0 0.0], 'k')
contour(D, d, funk_length, [0.01 0.01], 'k--')
% frequency
contour(D, d, funk_freq, [0.0 0.0], 'c')
contour(D, d, funk_freq, [0.01 0.01], 'c--')

%% objective & constaint function
function [spring_mass] = obj_function(smass)

spring_mass = smass;

end

%  constaint function
function [t1, f1, f2, L_spring, freq_s] = constraint(Tau2, F1, F2, L2, Lmin, freq1, Dv, p1, p2, nm)

%  calculating tau
tau12max = 600*10^6;
t1 = Tau2/tau12max - 1;

% calculating forces
% Nominal valve area:
Av = Dv^2*pi/4;
F1min = Av*p1;  % constant
F2min = Av*p2;  % constant
f1 = F1min/F1 -1;
f2 = F2min/F2 -1;

% calculating lengths, lmin is constant
L_spring = Lmin/L2 - 1;

% calculating freq
% nm is motor speed in rev per second, eigenfreqyency must not be equal to 
% a multiple of the groudn frequency so modulus(freq1/nm) / nm < 1
freq_s =  mod(freq1, nm) / nm - 1 ;
end


% Initialization
clf, hold off, clear
format long

% 1. Problem visualization
% Constant parameter values
springparams1;
w=1;
ktarget=10000; 
frtarget=300;

% Matrix of output values for combinations of design variables D and d 
D = [0.020:0.0005:0.040];
d = [0.002:0.00004:0.005];
for j=1:1:length(d)
  for i=1:1:length(D)
%   Analysis of valve spring.
    [svol,smass,bvol,matc,manc,Lmin,L2,k,F1,F2,Tau1,Tau2,freq1]=...
    springanalysis1(D(i),d(j),L0,L1,n,E,G,rho,Dv,h,p1,p2,nm,ncamfac,nne,matp,bldp);
 	 % Scaled objective function
     fobj(j,i) = ((k-ktarget)/ktarget)^2 + w*((freq1-frtarget)/frtarget)^2; 
     stiffness(j,i) = k;
     freq(j,i) = freq1;
  end
end

% Contour plot of scaled spring optimization problem
%contour(D, d, fobj,[0:0.05:0.2 0.2:0.1:0.5 0.5:0.5:2 2:5:100])
cc = [0.01 0.02 0.05];
contour(D, d, fobj,[cc 10*cc 100*cc 1000*cc 10000*cc 100000*cc 1000000*cc])
xlabel('Coil diameter D (m)'), ylabel('Wire diameter d (m)'), ...
   title('Figure 1     Spring stiffness and frequency optimization problem (w = 1.0)')
hold on
contour(D,d,stiffness,[10000 10000], 'Color', 'b')
contour(D,d,freq,[300 300], 'Color', 'k')
grid

% Initial design point:
xq = [0.022  0.0035];
plot(xq(1),xq(2),'o', 'Color', 'r');
hold on

% Initiation of optimization process:
again = 1;
cycle = 0;
D = xq(1);
d = xq(2);

options = optimset( 'Display', 'iter','TolX',1e-6, 'HessUpdate', 'bfgs', 'LargeScale', 'Off');
format long

[xnew,fval,exitflag] = fminunc(@(x) s_objw43(x, ktarget, frtarget, w), xq, options);

plot(xnew(1), xnew(2),'o', 'Color', 'r');
% Two variable valve spring problem - Exercise 3.1
% Visualization of SCALED spring stiffnes and frequency 
% optimization problem

% Initialization
clf, hold off, clear

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
     fobj(j,i) = abs((k-ktarget)/ktarget) + w*abs((freq1-frtarget)/frtarget);
     stiffness(j,i) = k;
     freq(j,i) = freq1;
  end
end
%% plotting
% Contour plot of scaled spring optimization problem
% contour(D, d, fobj,[0:0.05:0.2 0.2:0.1:0.5 0.5:0.5:2 2:5:100])
% cc = [0.01 0.02 0.05];
% contour(D, d, fobj,[cc 10*cc 100*cc 1000*cc 10000*cc 100000*cc 1000000*cc], 'showtext', 'on')
% xlabel('Coil diameter D (m)'), ylabel('Wire diameter d (m)'), ...
%    title('Figure 1: Spring stiffness and frequency optimization problem for w = 1.0')
% hold on
% contour(D,d,stiffness,[10000 10000],'showtext', 'off', 'Color','r')
% contour(D,d,freq,[300 300],'showtext', 'off', 'Color','k')
% grid




%% performing fminbnd
x1 = 0;
x2 = 10;

% start point and direction
xq = [0.022, 0.004];
sq1 = [0.002 0.0];
sq2 = [0.0 -0.0005];
sq3 = [0.002 -0.0005];
direction = sq2;

options = optimset('PlotFcns','optimplotfval', 'Display', 'iter','TolX',1e-8);
format long
[x,fval] = fminbnd(@(x) springobjw3(x, xq, direction, ktarget, frtarget, w), x1, x2, options);

% getting optimal point
xq2 = xq + x * direction;


%% getting arrows
% plot(xq2(1), xq2(2), 'r*') 
% quiver(xq(1), xq(2), sq1(1), sq1(2), 0)
% quiver(xq(1), xq(2), sq2(1), sq2(2), 0)
% quiver(xq(1), xq(2), sq3(1), sq3(2), 0)




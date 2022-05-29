% Valve spring design  - Exercise 5.1

% Computation of gradients of objective function and constraint g1.


% Initialization
clf, hold off, clear

% Note: Constant parameter values are read within the functions sprobj1 and sprcon1

% Design point for which gradients are computed 
% x = [0.024 0.004];
% x = [0.020522 0.003520];
% x = [0.02462 0.004035];
x = [0.022251 0.004075];


% Forward finite diffence gradients of objective function and constraints
hx = 1e-8; % vector of finite difference steps


% Objective function
fx = springobj1(x);
fx1plush = springobj1([x(1)+hx, x(2)]);
fx2plush = springobj1([x(1), x(2)+hx]);
dfdx1 = (fx1plush - fx)/hx;
dfdx2 = (fx2plush - fx)/hx;

% Constraints 
gx = springcon1(x);
gx1plush = springcon1([x(1)+hx, x(2)]);
gx2plush = springcon1([x(1), x(2)+hx]);
dgdx1 = (gx1plush - gx)./hx;
dgdx2 = (gx2plush - gx)./hx;

% setting up eq
div_g1 = [dgdx1(1), dgdx2(1)];
div_g2 = [dgdx1(2), dgdx2(2)];
div_g3 = [dgdx1(3), dgdx2(3)];
div_g4 = [dgdx1(4), dgdx2(4)];
div_g12 = cat(2, div_g1.', div_g2.');
div_g34 = cat(2, div_g3.', div_g4.');
div_g24 = cat(2, div_g2.', div_g4.');
div_g14 = cat(2, div_g1.', div_g4.');
div_f = [dfdx1;
        dfdx2];

mu34 = - inv(div_g34) * div_f
mu24 = - inv(div_g24) * div_f
mu12 = - inv(div_g12) * div_f
mu14 = - inv(div_g14) * div_f


 



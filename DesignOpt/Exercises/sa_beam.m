function sa_beam(n)
% sa_beam(n)
% n      = number of elements (optional, default 10)
%
% This script compares sensitivities computed by the Semi-Analytical
% approach for the beam example. Forward and central finite
% differences are considered. The output consists of two graphs with the
% relative error of the sensitivity results as a function of the relative
% design perturbation. Furthermore, the total evaluation time is listed.
% The number of elements n can be varied to examine the effect of mesh
% refinement.

% preparations
format long g



% define beam constants
L = 100;         % initial length
P = [1,1,0];     % load vector at tip: [Fx Fy M]
E = 2e5;         % Young's modulus
A = 2;           % cross-sectional area
I = 3;           % moment of inertia


% set number of elements
if ~nargin
    n = 10;  % relatively coarse mesh
    %n = 100; % fine mesh
end

% compute analytical sensitivity values, used as reference
dux_dL = P(1)/(E*A);
duy_dL = P(2)*L^2/(E*I);
du_ref = [dux_dL; duy_dL];

% relative design perturbations
pert = 10.^[-15:-1];

% compute forward SA sensitivities:
disp('Forward SA')
tic
sens_forward = fsa(pert, L, n, P, E, A, I);
toc

% compute central finite difference sensitivities:
disp('Central SA')
tic
sens_central = csa(pert, L, n, P, E, A, I);
toc

% plot the results
% ux sensitivity:
figure('color','w')
semilogx(pert, 100*(sens_forward(1,:)-dux_dL)/dux_dL,'ro--');
hold on
semilogx(pert, 100*(sens_central(1,:)-dux_dL)/dux_dL,'b*--');
% some additional lines
semilogx(pert, zeros(size(pert)),'k--');
semilogx(pert, ones(size(pert)),'k--');
semilogx(pert, -ones(size(pert)),'k--');

title(['Relative error of SA sensitivities for mesh with ',...
    num2str(n),' elements']);
xlabel('Relative design perturbation');
ylabel('Error [%] in du_x/dL')
legend('Forward SA','Central SA');

set(gca,'ylim',[-5 5]);

% uy sensitivity:
figure('color','w')
semilogx(pert, 100*(sens_forward(2,:)-duy_dL)/duy_dL,'ro--');
hold on
semilogx(pert, 100*(sens_central(2,:)-duy_dL)/duy_dL,'b*--');
% some additional lines
semilogx(pert, zeros(size(pert)),'k--');
semilogx(pert, ones(size(pert)),'k--');
semilogx(pert, -ones(size(pert)),'k--');

title(['Relative error of SA sensitivities for mesh with ',...
    num2str(n),' elements']);
xlabel('Relative design perturbation');
ylabel('Error [%] in du_y/dL')
legend('Forward SA','Central SA');

set(gca,'ylim',[-5 5]);



%% ----------------
% subroutines that perform the actual calculations

function sens = fsa (pert, L, n, P, E, A, I)
% ffd computes the forward semi-analytical sensitivity for every
% perturbation given in 'pert' and returns it in 'sens'.

% perform FEA for unperturbed case:
% (also returns nominal system matrix Knom)
[u,R,Kinv,Knom] = micro_fem(L,n,P,E,A,I);
% reorder displacement result into system vector:
Unom = u'; Unom = Unom(:);

% compute fsa derivative:
sens = zeros(2, length(pert));
for i=1:length(pert)
    Kpert = zeros(size(Knom));

    % apply forward difference to K matrix:
    % assemble perturbed K matrix:
    Lelem = (L*(1+pert(i))/n) *ones(1,n);
    for j=1:n
        % Calculate single element contribution
        Kelem = kelem_beam(Lelem(j), E, A, I);
        % Add it to the system matrix
        Kpert  = assemble_K(j,Kelem,Kpert);
    end    
    % compute pseudo-load (no design-dependent loads)
    pseudo_load = (-(Kpert - Knom)/(pert(i)*L))*Unom;
    % full sensitivity vector (apply boundary condition):
    allsens = Kinv*pseudo_load(4:end);
    
    % store only tip displacement sensitivity:
    sens(:,i) = [allsens(end-2); allsens(end-1)];
end


function sens = csa (pert, L, n, P, E, A, I)
% cfd computes the central semi-analytical sensitivity for every
% perturbation given in 'pert' and returns it in 'sens'.

% perform FEA for unperturbed case:
[u,R,Kinv,Knom] = micro_fem(L,n,P,E,A,I);
% reorder displacement result into system vector:
Unom = u'; Unom = Unom(:);

% compute fsa derivative:
sens = zeros(2, length(pert));
for i=1:length(pert)
    Kpert1 = zeros(size(Knom));
    Kpert2 = zeros(size(Knom));
    
    % apply central difference to K matrix:
    % assemble perturbed K matrix:
    Lelem = (L*(1-pert(i)))/n *ones(1,n);
    for j=1:n
        % Calculate single element contribution
        Kelem = kelem_beam(Lelem(j), E, A, I);
        % Add it to the system matrix
        Kpert1  = assemble_K(j,Kelem,Kpert1);
    end    
    Lelem = (L*(1+pert(i)))/n *ones(1,n);
    for j=1:n
        % Calculate single element contribution
        Kelem = kelem_beam(Lelem(j), E, A, I);
        % Add it to the system matrix
        Kpert2  = assemble_K(j,Kelem,Kpert2);
    end    
    
    % compute pseudo-load (no design-dependent loads)
    pseudo_load = (-(Kpert2 - Kpert1)/(2*pert(i)*L))*Unom;
    % full sensitivity vector (apply boundary condition):
    allsens = Kinv*pseudo_load(4:end);
    
    % store only tip displacement sensitivity:
    sens(:,i) = [allsens(end-2); allsens(end-1)];
end

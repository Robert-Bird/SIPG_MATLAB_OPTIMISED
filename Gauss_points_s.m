function [xsi,eta,w]=Gauss_points_s(p)

n=ceil((p+1)/2);
q=(1:n)./(1:2:2*n);
e=(1:n-1)./(3:2:2*n);
b0=2;
[~,D,V]=svd(diag(sqrt(q))+diag(sqrt(e),1));
x=(diag(D).^2)-1;
w=b0*V(1,:).^2;

gp_1=ones(length(x),1);
% face 1
xsi_f1 = -gp_1;
eta_f1 = -x;

% face 2
xsi_f2 = -x;
eta_f2 = gp_1;

% face 3
xsi_f3 = gp_1;
eta_f3 = x;

% face 4
xsi_f4 = x;
eta_f4 = -gp_1;

xsi=[xsi_f1;xsi_f2;xsi_f3;xsi_f4];
eta=[eta_f1;eta_f2;eta_f3;eta_f4];
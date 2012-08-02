% Find optimal H infinity controller
% coded A Sharma
%function [K, gamma_opt] = G2K(G);
function [K, gamma_opt] = G2K(G);

A = G.A; B = G.b; C = G.C; D = G.D;

% Calculate optimal controller

S = eye(size(D,2)) + D'*D;
R = eye(size(D,1)) + D*D';

[Y] = care((A-B/S*D'*C)', C', B/S*B', R);

[X] = care((A-B/S*D'*C), B, C'/R*C, S);

gamma_opt = sqrt(1 + max(eig(X*Y)));

F = S\(D'*C + B'*X);
W = (gamma_opt*gamma_opt - 1)*eye(size(X,2)) - Y*X;

Ac = W*(A-B*F) - gamma_opt*gamma_opt*(Y*C')*(C-D*F);
Bc = gamma_opt*gamma_opt*Y*C';
Cc = -B'*X;
Dc = -D';

% Get minimal realisation of controller
[U,EE,V] = svd(W);
o = rank(EE);
n = size(EE,1) - o;
E = EE(1:o,1:o);

Ac = U'*Ac*V;
Bc = U'*Bc;
Cc = Cc*V;
Dc = Dc;

% partition
a11 = Ac(1:o,1:o);
a12 = Ac(1:o,o+1:end);
a21 = Ac(1+o:end,1:o);
a22 = Ac(1+o:end,1+o:end);
b1 = Bc(1:o,:);
b2 = Bc(1+o:end,:);
c1 = Cc(:,1:o);
c2 = Cc(:,1+o:end);
d = Dc;

% introduce transformation & put minimal realisation
a = E\(a11-a12/a22*a21);
b = E\(b1-a12/a22*b2);
c = c1-c2/a22*a21;
d = d-c2/a22*b2;

K = ss(a,b,c,d);












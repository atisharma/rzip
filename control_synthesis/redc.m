% Reduce controller with closed-loop reduction method as
% described by El-Zobaidi, Jaimoukha, Limebeer.
% coded A Sharma
%
% expects weighted plant as input
%function [Kred, S, bound] = redc(G, order);
function [Kred, S, bound] = redc(G, order);

% Find suboptimal Hinf controller
grel = 1 + 1e-4;
[K, gamma, X, Y, Z, beta] = G2Ksubopt(G, grel);

% Find inverse weighted balancing transformation for theta
P = Z*Y*Z';
Q = X/(beta*beta*beta*beta);
%save test
[T,S] = XYbalance(P,Q,-1);

s = diag(S);

% Apply transformation
A = K.A; B = K.B; C = K.C; D = K.D;

A = T\A*T;
B = T\B;
C = C*T;
D = D;

% Truncate
A = A(1:order,1:order);
B = B(1:order,:);
C = C(:,1:order);
D = D;

m = 1 + 2.*s.*sqrt(1 + s.*s) + 2.*s.*s;
delta = -1 + prod(m(1+order:end));
bound = gamma/(1-delta*gamma);
if gamma*delta > 1
   bound = inf;
end

% return the reduced controller
Kred = ss(A,B,C,D);
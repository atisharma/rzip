% Find suboptimal H infinity controller
% coded A Sharma 2000
%function [K, gamma, X, Y, Z, beta] = G2Ksubopt(G, grel);
function [K, gamma, X, Y, Z, beta] = G2Ksubopt(G, grel);

A = G.A; B = G.b; C = G.C; D = G.D;

% Check D is zero
if (max(D)>0|min(D)<0) 
   error('Error: D is not zero');
   return
end

% Calculate optimal controller

[Y] = care(A', C', B*B');
[X] = care(A,  B,  C'*C);

gamma_opt = sqrt(1 + max(eig(X*Y)));
gamma = gamma_opt*grel;
beta = sqrt(1 - 1/(gamma*gamma));
Zinv = [eye(size(X,2)) - Y*X/(gamma*gamma*beta*beta)];
Z = inv(Zinv);

Ac = (A-B*B'*X) - Zinv\Y*C'*C/(beta*beta);
Bc = Zinv\Y*C';
Cc = -B'*X/(beta*beta);
Dc = -D';


K = ss(Ac,Bc,Cc,Dc);













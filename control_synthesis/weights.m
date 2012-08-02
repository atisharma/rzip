% Design weights for TCV plant
%
% function [Wo, Wi, Gw] = weights(G);
function [Wo, Wi, Gw] = weights(G);

% PS limits
%load ps

% OUTPUT WEIGHTS
%Dw = eye(size(G,1));
g0 = G.D - G.C/G.A*G.B;
[u,s,v] = svd(g0);
% generalised inverse
Dw = v*(s\u')*5e2;
Wo = ss(Dw);

% for shot 13333
% Have five outputs, weight ones with Ip by 1e-5.
%Dw = diag([2e5 .7e5 1e5 1 1]);
%Wo = ss(Dw);



% INPUT WEIGHTS
Wi = ss(eye(size(G,2)));
% or use inverse plant at DC
%Wi = ss(invplant(Wo*G));

% Add integrator, pole & zero
%w = ss(tf(1.5e7*[1 1e-1],[1 1e4 0]));
%Wi = w*Wi;

% Weighted plant is
Gw = Wo*G*Wi;

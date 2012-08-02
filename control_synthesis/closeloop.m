% Closes loop for G, K as described on p. 310, prob. 8.12
% Verified OK in simulink
function [H] = closeloop(G,K);

q = [G*K];

I  = ss(eye(size(q.D,1)));
Ig = ss(eye(size(G.D,1)));
Ik = ss(eye(size(K.D,2)));

J = inv([I - q]);
W = [G Ig];
E = [Ik; K];
H = E*J*W;

return

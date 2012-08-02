% Make input weighting matrix as generalised inverse of 
% output weighted plant

function [Wi] = invplant(h);

h0 = h.D + h.C/h.A*h.B;
% i.e. G = D + C*inv(sI - A)*B where s=0.
[u,s,v] = svd(h0);
Wi = tf(v/s*u');

return

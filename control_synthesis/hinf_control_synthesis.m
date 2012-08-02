% This script loads the tokamak model and p/s model,
% adds the weighting matrices to it for loop shaping,
% and synthesises a normalised left coprime uncertainty
% H-inf controller on the weighted plant.
% The closed loop transfer function is then found.
% Also the suboptimal controller is found
% and then reduced.
% Various figures may then be plotted.
% A Sharma 2000

gamma = inf;

[Wo,Wi,Gw] = weights(G);
% Gw = Wo*G*Wi;


[Kw,gamma] = G2K(Wo*G*Wi);
[Kwred, S] = redc(Gw, 17);
K = Wi*Kw*Wo;
Kred = Wi*Kwred*Wo;


%Hw = closeloop2(Gw,Kw);
%H = closeloop2(G,K);
%Hwred = closeloop2(Gw,Kwred);
%Hred = closeloop2(G,Kred);
%h = closeloopWdyn(G,K,Wi,Wo);

%He = esort(eig(H.A));
%%he = esort(eig(h.A));
%ge = esort(eig(G.A));
%e=esort(eig(Hw.A));
%er=esort(eig(Hwred.A));
%[e(1:5) er(1:5)]
%[norm(Hw,inf) norm(Hwred,inf)]

%makesim
%criteria
%criteria_red
%runsim

%save protcont Kred

return

%save controller



% Would like to add estimator (observer) matrices
% pulled from Alex Coutlis's / Parag Vyas's code 

% equilibrium according to Alex
R0 = 0.89;


load fls
load r_obs

A_r = [A_r zeros(1,18)];

Cmod = [[A_obs; A_r - R0*R0*A_obs(5,:)]];
Dmod = [[A_obs; A_r - R0*R0*A_obs(5,:)]];

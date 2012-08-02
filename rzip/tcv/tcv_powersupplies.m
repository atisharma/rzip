% TCV Power supplies as modelled by Jean-Yves Favez.
n_coils = max(coils);

% Power supply
taud=0.3e-3;
Aps=-eye(n_coils,n_coils)/taud;
Bps=eye(n_coils,n_coils)/taud;
Cps=eye(n_coils,n_coils);
Dps=zeros(n_coils,n_coils);

PSlimits = [1399 1399 648 648 648 648 648 648 648 648 1250 1250 1250 1250 1250 1250 1250 1250];
PScurrlimits = [31000 31000 7700 7700 7700 7700 7700 7700 7700 7700 7700 7700 7700 7700 7700 7700 7700 7700 ];

if max(coils) == 19;
% Time constant for g-coil
  tau_g = 0.2e-3;
  Aps(19,19) = -1/tau_g;
  Bps(19,19) = 1/tau_g;
  PSlimits = [PSlimits 566];
  PScurrlimits = [PScurrlimits 2e3];
end 

Gps = ss(Aps,Bps,Cps,Dps);


save ps Gps PSlimits PScurrlimits
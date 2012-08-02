% Set up simulation

delta_t = 1e-4;

% Make discrete time simulink model
Gd	= c2d(G,delta_t,'tustin');
Gwd	= c2d(Gw,delta_t,'tustin');
Kwd	= c2d(Kw,delta_t,'tustin');
Kwredd	= c2d(Kwred,delta_t,'tustin');

%Hd	= c2d(H,delta_t,'foh');
%hd	= c2d(h,delta_t,'foh');

Wid	= c2d(Wi,delta_t,'tustin');
Wod	= c2d(Wo,delta_t,'tustin');


PSlimits = [1399 1399 648 648 648 648 648 648 648 648 1250 1250 1250 1250 1250 1250 1250 1250 566];
PSlimits = PSlimits(1:size(K,1));

% Make the signals [p_vert tri_out tri_in zIp Ip Psi_R]
T = [0:delta_t:.6]';
ref = zeros(length(T),5);
ref(find(T>0.35 & T<0.4),1) = .02;
ref(find(T>0.15 & T<0.2),2) = .02;
ref(find(T>0.25 & T<0.3),3) = .02;
ref(find(T>0.05 & T<0.1),4) = 2000;
ref(find(T>0.45 & T<0.5),5) = 5000;
%ref(1350:1500,5) = 2e3;


save simmat

%% Closes loop for G, K with STATIC weighting Wi, Wo.
%% Inputs  [[r' w' v'] u']'
%% Outputs [[y' u'] e']'
%% So controller closes u to e.
%% r = ref
%% w = dist. on inputs
%% v = dist. on outputs
%% u = controller output
%% y = plant output + disturbance
%% e = plant output - ref.

function [H] = closeloopw(G,K,Wi,Wo);

A = G.A; B = G.B; C = G.C; D = G.D;
a = K.A; b = K.B; c = K.C; d = K.D;
Di = Wi.d; Do = Wo.d;

x = size(A,1);
y = size(D,1);
u = size(D,2);

save cli

B1 = [B*Di zeros(x,y)];
B2 = B;
C1 = [Do*C; zeros(u,x)];
C2 = C;
D11 = [Do*D*Di eye(y,y); zeros(u,u) zeros(u,y)];
D12 = [Do*D; eye(u,u)];
D21 = [D*Di inv(Do)];
D22 = D;


M = inv(eye(size(D22*d)) - D22*d);

Ah = [(A + B2*d*M*C2) (B2*(eye(size(d*M*D22)) + d*M*D22)*c); b*M*C2 (a+b*M*D22*c)];
Bh = [B1 + B2*d*M*D21; b*M*D21];
Ch = [(C1 + D12*d*M*C2) (D12*(eye(size(d*M*D22))+ d*M*D22)*c)];
Dh = [D11 + D12*d*M*D21];


H = ss(Ah, Bh, Ch, Dh);


return

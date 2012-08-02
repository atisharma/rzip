% Closes loop for G, K as described on p. 137

function [H] = closeloop2(G,K);

A = G.A; B = G.B; C = G.C; D = G.D;
a = K.A; b = K.B; c = K.C; d = K.D;

x = size(A,1);
y = size(D,1);
u = size(D,2);

save cli
% Form P
B1 = [B zeros(x,y)];
B2 = B;
C1 = [zeros(u,x); C];
C2 = C;
D11 = [zeros(u,u) zeros(u,y); D eye(y,y)];
D12 = [eye(u,u); D];
D21 = [D eye(y,y)];
D22 = D;


M = inv(eye(size(D22*d)) - D22*d);

Ah = [(A + B2*d*M*C2) (B2*(eye(size(d*M*D22)) + d*M*D22)*c); b*M*C2 (a+b*M*D22*c)];
Bh = [B1 + B2*d*M*D21; b*M*D21];
Ch = [(C1 + D12*d*M*C2) (D12*(eye(size(d*M*D22))+ d*M*D22)*c)];
Dh = [D11 + D12*d*M*D21];


H = ss(Ah, Bh, Ch, Dh);


return

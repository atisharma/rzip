%% Closes loop for G, K with dynamic weighting Wi, Wo.
%% Inputs  [[w' v'] u_k']'
%% Outputs [[y' u'] y_k']'
%% So controller closes u to y.
%% w = dist. on inputs   --- 
%% v = dist. on outputs  --- All these are weighted
%% u = controller output ---
%% y = plant output + disturbance ---
%% 
%% Controller K acts on y_k to give u_k
%% K = Wi*Kw*Wo
%% Gw = Wo*G*Wi
%% Kw designed for Gw

function [H] = closeloopWdyn(G,K,Wi,Wo);

t = Wo*G*Wi;
x = size(t.A,1);
y = size(t.D,1);
u = size(t.D,2);

% Put in general form (P)

Iu = ss(eye(u));
Iy = ss(eye(y));

P = [Wo*G*Wi        Iy                   Wo*G; 
     ss(zeros(u,u)) ss(zeros(u,y))       Iu;
     G*Wi           inv(Wo)*Iy           G];

A = P.A; B = P.B; C = P.C; D = P.D;
a = K.A; b = K.B; c = K.C; d = K.D;

% Partition the P matrix

B1 = B(:,           1:u+y);
B2 = B(:,           u+y+1:end);
C1 = C(1:y+u,       :);
C2 = C(y+u+1:end,   :);
D11 = D(1:y+u,      1:y+u);
D12 = D(1:y+u,      y+u+1:end);
D21 = D(y+u+1:end,  1:y+u);
D22 = D(y+u+1:end,  1+y+u:end);

% Give minimal realisation of lower LFT

M = inv(eye(size(D22*d)) - D22*d);

Ah = [(A + B2*d*M*C2) (B2*(eye(size(d*M*D22)) + d*M*D22)*c); b*M*C2 (a+b*M*D22*c)];
Bh = [B1 + B2*d*M*D21; b*M*D21];
Ch = [(C1 + D12*d*M*C2) (D12*(eye(size(d*M*D22))+ d*M*D22)*c)];
Dh = [D11 + D12*d*M*D21];


H = ss(Ah, Bh, Ch, Dh);
save clint

return

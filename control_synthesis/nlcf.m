% Find normalised left coprime factorisation of G of (ABCD).
% A Sharma 8/11/00
% function [N, M] = nlcf(A,B,C,D);

function [N, M] = nlcf(A,B,C,D);


% We first solve the Ricatti equation described on P.420 of [1]
% Assume G is stabilisable, detectable.
S = D'*D+eye(size(D,2));
R = D*D'+eye(size(D,1));
Y = care([A-(B/S)*D'*C]', C', [(B/S)*B'], R);
H = (B*D' + Y*C')/R;
F = R^-0.5;

sys = ss( [A-H*C], [B-H*D H], [-F*C], [F*D F] );
whos
[Q] = tf(sys);

% Q = [N M], size n-by-(inputs+outputs)
N = Q(:,size(G,2));
M = Q(:,(size(G,2)+1):end);
% Check
size(M)
size(N)

return

%------------------------------------------------------


% [1] Green & Limebeer, 'Linear Robust Control'
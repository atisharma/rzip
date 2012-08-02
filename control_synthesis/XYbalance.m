% Uses code from Imad Jaimoukha
function [T,S]=XYBalance(X,Y,p)
%        [T,S]=XYBalance(X,Y,p)
%   Ti*X*Ti'=T'*Y*T=S=diag(diag(S))
%    p=1: increasing magnitude
%    p=-1: decreasing magnitude

Xc = chol(X)';
XcYXc = Xc'*Y;
x = [];
n = size(X,1);

for i = 1:n, x(1:i,i) = XcYXc(1:i,:)*Xc(:,i); end;
XcYXc = triu(x) + triu(x,1)';
[V,S] = schur(XcYXc);
if p == -1, [x,i] = sort(-diag(S)); S = S(i,i); V = V(:,i);
elseif p == 1, [x,i] = sort(diag(S)); S = S(i,i); V = V(:,i);
end

s = diag(S);
s(find(s<eps)) = eps;
S = diag(s);
S = diag(sqrt(abs(diag(S))));
T = Xc*V*diag(diag(S).^(-.5));

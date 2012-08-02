%% Load TCV model, use Schur model reduction to remove 3 zero eigenvalues
% see IEEE Transactions on Control Systems Technology Vol. 13 No. 3 (2005)
% for the significance of this
% A S Sharma 2000
function [ar,br,cr,dr] = schurred(A,B,C,D);

% There are typically 3 zero eigenvalues for a tokamak
n=size(A,1);nn=n-rank(A);nsa=n-nn;

% Partition realisation of system such that last three states are associated with these modes
[V,a]=blkrsch(A,5,nsa);
A11=a(1:nsa,1:nsa);A12=a(1:nsa,nsa+1:n);A22=a(nsa+1:n,nsa+1:n);
A21=a(nsa+1:n,1:nsa);
b=V'*B;c=C*V;
B1=b(1:nsa,:);B2=b(nsa+1:n,:);
C1=C(:,1:nsa);C2=c(:,nsa+1:n);
d=D;
% D is zeros usually

% A21, A22 near zero but not that near zero.
X=lyap(A11,-A22,A12);

ar=[A11 (A11*X - X*A22 +A12); zeros(nn,nsa) A22];
T = [eye(size(X,1)) -X;zeros(size(X')) eye(size(X,2))];
Td = [eye(size(X,1)) X;zeros(size(X')) eye(size(X,2))];
ar = T*a*Td;
br = T*b;
cr = c*Td;
dr=D;

% Now truncate the three last states
ar = ar(1:nsa,1:nsa); br = br(1:nsa,:); cr = cr(:,1:nsa);
% Could consider leaving the poles in and doing pole shifting to add
% integral action to controller


% Check the sv's are the same

%[sv, w] = sigma(A,B,C,D);
%svs = sigma(a,b,c,d);
%[svr] = sigma(ar,br,cr,dr,w);
%svsz = sigma([A11 A12; zeros(size(A21)) A22], b,c,d);
%loglog(w,svs,'-', w,svr,'o', w,svsz,'+')


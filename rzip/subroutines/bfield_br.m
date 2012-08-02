%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to define the radial b field due to a circular filament
% a is radius of filament, b is radius of position, 
% c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function br = bfield_br(a,b,c,cur)

% set up constants
 mu0=4e-07*pi;

% ensure matrices are the same orientation
 if(size(a,1)  <size(a,2))  ; a=a'; end;
 if(size(b,1)  <size(b,2))  ; b=b'; end;
 if(size(c,1)  <size(c,2))  ; c=c'; end;
 if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
 com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;

% construct the elliptic integrals
 k2 =  4*a.*b./((a+b).^2 + c.^2);
 [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;

 f=com*1./sqrt((a+b).^2+c.^2);
 f1=-K+(a.*a+b.*b+c.*c)./((a-b).^2.+c.^2.).*E;

% define the value of br (Smythe pg 290)
 br=2.e-7*cur.*c.*f./b.*f1;

 clear E K f f1 k2 mu0;

 return;

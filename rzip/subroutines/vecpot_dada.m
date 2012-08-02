%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to define the radial gradient of the vector potential
% due to a movement of a circular filament
% a is radius of filament, b is radius of position, 
% c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dAda = vecpot_dAda(a,b,c,cur)

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
 k = sqrt(k2);
[K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;
 dEdk= (E./k-K./k);
 dKdk= (E./((1-k.^2).*k)-K./k);

 r=b;
 z=c;
 ak=k;

    g3=(a+r).^2+z.^2;
 dg3da=2*(a+r);

 dkda=(4*r.*g3-4*dg3da.*a.*r)./(2*k.*g3.^2);


 F3=(2-k2).*K-2*E;
 dF3da=(2-k2).*dKdk.*dkda - 2*k.*dkda.*K - 2*dEdk.*dkda;

 dAda= 2e-7*cur./sqrt(r).*( sqrt(a).*dF3da./k + F3./(2*sqrt(a).*k) - sqrt(a).*dkda.*F3./k2 );

 clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to define the radial gradient of the radial magnetic field
% due to a movement of a circular filament
% a is radius of filament, b is radius of position, 
% c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbrda = bfield_dbrda(a,b,c,cur)

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

    g2=(a-r).^2+z.^2;
 dg2da=2*(a-r);

    g3=(a+r).^2+z.^2;
 dg3da=2*(a+r);

    g4=a.^2+r.^2+z.^2;
 dg4da=2*a;


 dkda=(4*r.*g3-4*dg3da.*a.*r)./(2*k.*g3.^2);


 F3=-K+g4.*E./g2;
 F4=g3.^(-1/2);

 dF3da=-dKdk.*dkda+(dg4da.*E+dEdk.*dkda.*g4)./g2 - dg2da.*g4.*E./(g2.^2);
 dF4da=-dg3da./(2*g3.^(3/2));

 dbrda= 2e-7*cur.*z./r.*(dF3da.*F4+dF4da.*F3);

 clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

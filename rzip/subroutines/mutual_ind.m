%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file to compute the mutual inductance between two coaxial circular filaments
% input parameters 	a,b	=	radius of coils
%			c	=	separation of planes of coils
% ref. Smythe p.335
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function mut_ind = mutual_ind(a,b,c)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
 if(size(a,1)  <size(a,2))  ; a=a'; end;
 if(size(b,1)  <size(b,2))  ; b=b'; end;
 if(size(c,1)  <size(c,2))  ; c=c'; end;

% set up a dummy mapping array
 com=ones(size(a));


% construct the elliptic integrals
 k2 =  4*a.*b./((a+b).^2 + c.^2);
 k = sqrt(k2);
 i_ch=find(k2==1);
 k2(i_ch)=0.1*ones(size(i_ch));
[K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;

 mut_ind = 2*mu0*com./k.*sqrt(a.*b).*((com*1.-0.5*k2).*K - E);
 mut_ind(i_ch)=zeros(size(i_ch));

% mut_ind(find(isnan(mut_ind)))=zeros(size(find(isnan(mut_ind))));
% mut_ind(find(isinf(mut_ind)))=zeros(size(find(isinf(mut_ind))));

 return;


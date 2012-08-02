%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the self inductance of a coil of square
% cross section, width w height h and major radius a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [sel_ind] = self_sqre(a,w,h)

 if(size(a,1)  <size(a,2))  ; a=a'; end;
 if(size(w,1)  <size(w,2))  ; w=w'; end;
 if(size(h,1)  <size(h,2))  ; h=h'; end;

 mu0=4e-07*pi;

 sel_ind = mu0*a.*(log(8*a./(w+h))-.5);

 clear aphi mu0 k;

 return;

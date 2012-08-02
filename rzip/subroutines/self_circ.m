%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the self inductance of a coil of circular
% cross section, radius r and major radius a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [sel_ind] = self_circ(a,r)

 if(size(a,1)  <size(a,2))  ; a=a'; end;
 if(size(r,1)  <size(r,2))  ; r=r'; end;

 mu0=4e-07*pi;

 sel_ind = mu0*a.*(log(8*a./r)-1.75);

 clear aphi mu0 k;

 return;

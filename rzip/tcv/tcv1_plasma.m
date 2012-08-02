% a function to create equilibrium data for a Machine

function [J_PL,PS_PL,A_PL,R_PL,Z_PL,BETA_P,LI,I_AC]=TCV1_PLASMA(shot,time);

    
 eval(['load plasma_eqm_tcv_'num2str(shot),'.mat;']);
 LI=1.3;

 load /home/ps/sunfs1/d4/ati/matlab/model/psi13333
 PS_PL=PSI(:);

 J_PL = J_PL(:);
 A_PL = ones(size(J_PL));
 DA   = abs((R_PL(2)-R_PL(1))*(Z_PL(2)-Z_PL(1)));
 A_PL = A_PL*DA;

 NR=length(R_PL);NZ=length(Z_PL);N_grid_points=NR*NZ;
 RG=zeros(N_grid_points,1);ZG=zeros(N_grid_points,1);JG=zeros(N_grid_points,1);
 for k=1:NZ;
  ZG((k-1)*NR+1:(k-1)*NR+NR)=Z_PL(k)*ones(NR,1);
 end; 
 for k=1:NR;RG(k:NR:N_grid_points)=R_PL(k)*ones(NZ,1); end;
 JG=J_PL(:);Z_PL=ZG;R_PL=RG;

return;

% a function to take an array of RZ positions and interpolate two
% intermediate points between each position based on a parabola 
% throught these points and the one after. For a closed line
% where the first and last point are the same.

function [RR,ZZ,XX,YY,NN]=add_points_closed(fil_nam,NP);

 eval(['load ' fil_nam]);

% define the number of points 
 NN=max(size(RC)); 
 RCT=[RC RC(2:5)];
 ZCT=[ZC ZC(2:5)];

% define a large array to keep track of points 
 BX1=[];BY1=[];

% begin loop

 for k=2:NN;

% first add current point to B arrays
  BX1=[BX1 RC(k-1)];BY1=[BY1 ZC(k-1)];

% define indices for this and next two points
  k0=k-1;k1=k;k2=k+1;

% determine the local gradient. If too steep need to consider Z as variable
    if((RCT(k1)-RCT(k0))); 
   grad_line=(ZCT(k1)-ZCT(k0))/(RCT(k1)-RCT(k0));
  else; grad_line=1000; end;
  if(abs(grad_line)<1)  
   X_test=[RCT(k0) RCT(k1) RCT(k2)];
   Y_test=[ZCT(k0) ZCT(k1) ZCT(k2)];
   V_mark=1;
  else
   X_test=[ZCT(k0) ZCT(k1) ZCT(k2)];
   Y_test=[RCT(k0) RCT(k1) RCT(k2)];
   V_mark=0;
  end;

% fit parabola to three points
   M=[X_test(1)^2 X_test(1) 1;
      X_test(2)^2 X_test(2) 1;
      X_test(3)^2 X_test(3) 1];
   coeffs=inv(M)*[Y_test(1);Y_test(2);Y_test(3)];

% define intermediate points
   dx=(X_test(2)-X_test(1))/(NP+1); 
   XIN=X_test(1)+[1:NP]*dx;
   YIN=coeffs(1)*XIN.^2+coeffs(2)*XIN+coeffs(3);
 
% according to local gradient XIN and YIN could be R or Z 
 if(V_mark);
   BX1=[BX1 XIN];BY1=[BY1 YIN];
  else
   BX1=[BX1 YIN];BY1=[BY1 XIN];
  end;

 end; % finish loop

% finish off array with start point
   BX1=[BX1 RC(1)];BY1=[BY1 ZC(1)];

   RR=BX1; ZZ=BY1; 
   NN=max(size(RR))-1;

   XX=0.5*(RR(1:NN)+RR(2:NN+1));  YY=0.5*(ZZ(1:NN)+ZZ(2:NN+1));
   NN=max(size(XX)); 

 clear RCT ZCT RC ZC BX1 BY1 k M coeffs dx XIN YIN X_test Y_test V_mark;
 clear k0 k1 k2 grad_line;

 return;

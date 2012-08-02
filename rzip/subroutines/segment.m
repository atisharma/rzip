% given a rectangular conductor segment into NP smaller conductors along
% its longest side. R,Z = centre coordinates, DR,DZ = half width.

 function [Ro,Zo,DRo,DZo] = segment(R,Z,DR,DZ,NP);

% find longest side
  Rlong=1;if(DZ>DR);Rlong=0;end;

% use X and Y as intermediate, X is discretised

  if(Rlong)
   X = R; DX = DR; Y = Z; DY = DZ;
  else
   X = Z; DX = DZ; Y = R; DY = DR;
  end;

% define limits
  Xwid = 2*DX/NP;
    X0 = X-DX+Xwid/2;
  
  if(Rlong)
   Ro = X0+[0:NP-1]'*Xwid ; DRo = Xwid*ones(NP,1); 
   Zo = Y*ones(NP,1);       DZo = 2*DY*ones(NP,1);
  else
   Zo = X0+[0:NP-1]'*Xwid ; DZo = Xwid*ones(NP,1); 
   Ro = Y*ones(NP,1);       DRo = 2*DY*ones(NP,1);
  end

return

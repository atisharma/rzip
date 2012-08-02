%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% link filaments to coils in the appropriate manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [new_data]= coil_from_filaments(old_data,link,nc,nf1,nf2);

% data is in order: coils then VV. nc is number
% of coils. nf1 is number of filaments for the coils
% nf2 is the number of filaments for the VV
% link  is the coil linkage information
  
 new_data=[]; 

% for the active coils
 for k=1:nc
  if(length(link(k,1):link(k,2))~=1)
   new_data    =[new_data; sum(old_data(link(k,1):link(k,2),:))];
  else
   new_data    =[new_data; old_data(link(k,1):link(k,2),:)];
  end;
 end;

% add the passive data
 new_data   =[new_data;old_data(nf1+1:nf1+nf2,:)];

 return;


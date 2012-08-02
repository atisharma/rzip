%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a function to ensure vectors are columns;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [A]=check_is_column(A);

% only correct vectors
 if(any(size(A)==1))
  if(size(A,1)<size(A,2));A=A';end;
 end;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       unique_number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a file which give a column vector of numbers, counts the different numbers and defines the
% range in the list of each name (index range) and total length of vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [nn,nf,index_range,names]=unique_number(numbers);

% define total length and create a row of the names
  nf=size(numbers,1); 

  numbers=sort(numbers);

% count the number of different names 
  nn=1; test_num=numbers(1);
	
  names(nn)=test_num;
  ind=find(numbers==test_num);
  index_range=[min(ind) max(ind)];
  while(max(ind)~=nf);
   test_num=numbers(max(ind)+1);
   ind=find(numbers==test_num);
   names=[names;test_num];
   index_range=[index_range;[min(ind) max(ind)]]; 
   nn=nn+1;
  end;

return;



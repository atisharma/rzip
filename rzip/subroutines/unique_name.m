% a file which give a colum vector of text names
% counts the different names(nn) and defines the
% range in the list of each name (index range)
% and total length of string

 function [nn,nf,index_range,names]=unique_name(name_string);

% define total length and create a row of the names
  nf=size(name_string,1); 
  nc=size(name_string,2); 

% a potential error is if a name is made of all the same characters!
% To avoid this scenario add a dummy separation variable
  name_string2=[name_string,ones(nf,1)*'*'];
  name=name_string2';name=name(:)';

% count the number of different names 
 nn=1; test_name=name_string(1,1:nc);
 
 names(nn,1:nc)=test_name;
 ind=[findstr(name,test_name)+nc]/(nc+1);
 index_range=[min(ind) max(ind)];
 while(max(ind)~=nf);
  test_name=name_string(max(ind)+1,1:nc);
  ind=[findstr(name,test_name)+nc]/(nc+1);
  names=[names;test_name];
  index_range=[index_range;[min(ind) max(ind)]]; 
  nn=nn+1;
 end;

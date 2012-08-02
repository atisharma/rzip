% a script to print out the lines relevant to getting PROTEUS
% working for an equilibrium

  function for_pro(shot)

% first load in the appropriate data
cd data
  eval(['load shot_',num2str(shot)]);

% load in the plasma mesh
[nel,nno,link,code,noder,nodez,ni]=readmesh('C:/PROTEUS/TCV/TCV_LST/fort.8');
% see which one?

% evaluate the mutual coupling from all points with current
% to all grid points!

% all the conductors
  code_els = [ones(16,1)*1;ones(6,1)*2;[3:18]';ones(6,1)*19];
  code_fid = [ones(16,1)*1;ones(6,1)*104/54;ones(8,1)*34/9;ones(8,1)*4;ones(6,1)*6/54];
  coil_node =[];
  coil_ind1 =[];
  coil_ind2 =[];
  coil_fidd =[];
  for k=1:44
   find_coil = find(code == k);
   link_nodes = link(find_coil,:);link_nodes=link_nodes(:); 
   [nn,nf,index_range,nums]=unique_number(link_nodes);
   coil_node = [coil_node; nums];
   coil_ind1 = [coil_ind1; ones(size(nums))*k];
   coil_ind2 = [coil_ind2; ones(size(nums))*code_els(k)];
   coil_fidd = [coil_fidd; ones(size(nums))*code_fid(k)];
  end;


% find plasma nodes within 2 cm of the plasma rzip nodes

  find_plasma = find(code == -1);
  plasma_nodes = [];
  for k=1:length(find_plasma)
   plasma_nodes = [plasma_nodes ; link(find_plasma(k),:)];
  end;
  [Nun, d1,d2,plasma_nodes] = unique_number(plasma_nodes(:));

  valid_node =[];
  valid_curr =[];
  fJ = find(J_PL);
  for k=1:Nun
   dist  = sqrt( (R_PL(fJ)-noder(plasma_nodes(k))).^2 + (Z_PL(fJ)-nodez(plasma_nodes(k))).^2);
   if(any(dist<.02));
    [dummy,idist] = min(dist); 
    valid_node = [valid_node;k]; 
    valid_curr = [valid_curr;J_PL(fJ(idist))*A_PL(fJ(idist))];
   k=k+1; end;
  end;
%  plot(noder(plasma_nodes),nodez(plasma_nodes),'b+',noder(plasma_nodes(valid_node)),nodez(plasma_nodes(valid_node)),'ro')

% so the total nodes with a current are 
  total_nodes = [coil_node;plasma_nodes(valid_node)];
  N_tn = length(total_nodes);
% and those currents are
  total_currs = zeros(N_tn,1);
  for k = 1:length(coil_node)
   total_currs(k) = I_AC(coil_ind2(k))*coil_fidd(k);
  end;
  total_currs([length(coil_node)+1]:N_tn) = valid_curr;


  bound = find(ni==2 | ni==1 | ni==3);
  non_bound = find(ni~=2 & ni~=1 & ni~=3)';
  LNB = length(non_bound);

% need to map these to whole grid
  L = zeros(LNB,N_tn);
  for k=1:N_tn
   L(:,k)=mutual_ind(noder(total_nodes(k)),noder(non_bound),nodez(non_bound)-nodez(total_nodes(k)));
  end;

  for k=1:N_tn
   L(k,k)=self_sqre(noder(total_nodes(k)),.02,.02);
  end;

  PS = zeros(size(noder));

  PS(non_bound) = L*total_currs;
  
  
  fid = fopen(['shot_psi_',num2str(shot)],'w');
  for k=1:length(noder)
   fprintf(fid,'%13.8e\n',PS(k));
  end;
  fprintf(fid,'%13.8e\n',max(PS));
  fprintf(fid,'%13.8e\n',PSI_0);
  fclose(fid);
  
  cd ..
  
  return;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                 mut_ind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file to compute the mutual inductance between two coaxial circular filaments
% input parameters 	a,b	=	radius of coils
%			c	=	separation of planes of coils
% ref. Smythe p.335
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% self_sqre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the self inductance of a coil of square cross section, width w, height h and radius a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [sel_ind] = self_sqre(a,w,h)

  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(w,1)  <size(w,2))  ; w=w'; end;
  if(size(h,1)  <size(h,2))  ; h=h'; end;

  mu0=4e-07*pi;

  sel_ind = mu0*a.*(log(8*a./(w+h))-.5);

  clear aphi mu0 k;

 return;




  PS_EXT = L*[total_currs(1:length(coil_node));zeros(size(valid_curr))];
  plot3(noder,nodez,PS,'ro')

  load /home/ps/sunfs1/d4/jw1/Matlab/rzip/tcv/tcv_plasma_full_v3
  PSI_PL_EXT = mut_struct_grid(1:19,:)'*I_AC(1:19);


  plot3(noder(plasma_nodes),nodez(plasma_nodes),PS(plasma_nodes),'ro')
  plot3(noder(plasma_nodes(valid_node)),nodez(plasma_nodes(valid_node)),valid_curr,'ro')





  
% symmetric amps
  CPS = [3 10;4 9;5 8;6 7;11 18;12 17;13 16; 14 15];
  amp =I_AC;
  for k=1:8
   a_sym = mean(I_AC(CPS(k,:)));
   amp(CPS(k,1))=a_sym;
   amp(CPS(k,2))=a_sym;
  end


  disp(['   ',num2str(max(max(PSI_PL+PSI_0))/2/pi)])
  disp(['   ',num2str(PSI_0/2/pi)])

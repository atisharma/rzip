% TCV description with conducting tiles, 

function [R_AC,Z_AC,RES_AC,DR_AC,DZ_AC,LA_AC,NT_AC,R_PA,Z_PA,RES_PA,DR_PA,DZ_PA,LA_PA,R_BP,Z_BP,TH_BP,LA_BP,R_FL,Z_FL,LA_FL]=tcv_tiles();

% Set up the maximum number of eigenmodes used array in here
  global N_eigen_modes_array_max

  N_eigen_modes_array_max = [40 40];


% load ./tcv_magnetics;       %(original file sent for distribution)
% load ./Basic_data;          %(has no filaments or cross sections, this does).
  load tcv_structure;         %(relevant numbers extracted from above)



% ACTIVE COILS
%%%%%%%%%%%%%%%%

  used_fils=1:830;
  g_coils=825:830;
  no_co=29;

 R_AC=fil_r(used_fils);
 Z_AC=fil_z(used_fils);
 RES_AC=2*pi*fil_re*R_AC(used_fils)./fil_cs(used_fils);
 DR_AC=fil_dr(used_fils);
 DZ_AC=fil_dz(used_fils);

% no DR and DZ for gcoils so bluff them
  DR_AC(g_coils)=sqrt(fil_cs(g_coils)/pi);
  DZ_AC(g_coils)=sqrt(fil_cs(g_coils)/pi);


% add the missing coils to nt_c (could have used nt_a*nt_b)
% but necessary to go from filament list to coils
 nt_adj=[1 3 3 0 0 0 0 2 2 2 2 2 2 2 2]';
 nt_c(1:15)=nt_c(1:15)+nt_adj;

% define an array which matches each coil to a range of filaments
 nt_i(1,1)=1;nt_i(1,2)=nt_c(1);
 for k=2:no_co;
  nt_i(k,1)=1+sum(nt_c(1:k-1));
  nt_i(k,2)=sum(nt_c(1:k));
 end;


% can no define the label array
 LA_AC=[]; 
% OH1
 for k=1:1
  ffk='00';ffl='00';
  for l=1:nt_c(1);
   if(k>=10&k<99);ffk='0';end;  if(k>=100);ffk='';end; 
   if(l>=10&l<99);ffl='0';end;  if(l>=100);ffl='';end; 
%   eval(['LA_AC=[LA_AC;''A001F',ffl,num2str(l,3),'''];']);
   eval(['LA_AC=[LA_AC;''OH101F',ffl,num2str(l,3),'''];']);
  end;
 end;
% OH2
 for k=2:7
  ffk='00';ffl='00';
  for l=1:nt_c(k);
   if(k>=10&k<99);ffk='0';end;  if(k>=100);ffk='';end; 
   if(l>=10&l<99);ffl='0';end;  if(l>=100);ffl='';end; 
%   eval(['LA_AC=[LA_AC;''A002F',ffl,num2str(l,3),'''];']);
   eval(['LA_AC=[LA_AC;''OH20',num2str(k-1),'F',ffl,num2str(l,3),'''];']);
  end;
 end;
 % for E  
 for k=8:15;
  ffk='00';ffl='00';
  for l=1:nt_c(k);
   ki=k-5;
   if(ki>=10&ki<99);ffk='0';end;  if(ki>=100);ffk='';end;
   if(l>=10&l<99);ffl='0';end;  if(l>=100);ffl='';end; 
%   eval(['LA_AC=[LA_AC;''A',ffk ,num2str(ki,3),'F',ffl,num2str(l,3),'''];']);
   eval(['LA_AC=[LA_AC;''E_',num2str(k-7),'01F',ffl,num2str(l,3),'''];']);
  end;
 end;
 % for F
 for k=16:23;
  ffk='00';ffl='00';
  for l=1:nt_c(k);
   ki=k-5;
   if(ki>=10&ki<99);ffk='0';end;  if(ki>=100);ffk='';end;
   if(l>=10&l<99);ffl='0';end;  if(l>=100);ffl='';end; 
%   eval(['LA_AC=[LA_AC;''A',ffk ,num2str(ki,3),'F',ffl,num2str(l,3),'''];']);
   eval(['LA_AC=[LA_AC;''F_',num2str(k-15),'01F',ffl,num2str(l,3),'''];']);
  end;
 end;
% the gcoils
 ffk='00';ffl='00';
 for k=1:6;
%  eval(['LA_AC=[LA_AC;''A019F',ffl,num2str(k,3),'''];']);
  eval(['LA_AC=[LA_AC;''G_10',num2str(k),'F00',num2str(k,3),'''];']);
 end;

% TCV FIDDLES
% make the TCV 'corner corrections'

 NT_AC=[ones(824,1)];
 NT_AC=[NT_AC(1:144,:)*143/144;
        NT_AC(145:208,:)*29/32;
        NT_AC(209:248,:);
        NT_AC(249:536,:)*34/36;
        NT_AC(537:824,:);
        ones(3,1);
       -ones(3,1)];

% then adjust resistivity  to get the res_a values! hmmmm
 RES_AC(1:144)  = (9.930567081532160e-01)^2*RES_AC(1:144);
 RES_AC(145:248)= (9.534982912147054e-01)^2*RES_AC(145:248);
 RES_AC(249:536)= (9.444447551840814e-01)^2*RES_AC(249:536);
 
% correct gcoils again to get the self inductance right
 ff=.288;
 DR_AC(g_coils)=DR_AC(g_coils)/ff;
 DZ_AC(g_coils)=DZ_AC(g_coils)/ff;


% PASSIVE STRUCTURE
 R_PA=vv_r;
 Z_PA=vv_z;
 RES_PA=res_v;
 DR_PA=ves_w;
 DZ_PA=ves_h;
 LA_PA=vvs;

 LA_PA=[];ff='00';
 for k=1:length(R_PA);
  if(k>=10&k<99);ff='0';end;  if(k>=100);ff='';end;
  eval(['LA_PA=[LA_PA;''VV',ff,num2str(k,3),'F001''];']);
 end;

% add the tiles

 load /home/ps/sunfs1/d4/jw1/Matlab/rzip/tcv-tiles/tiles
%   tile_range=[1:21 34:61];
%   tile_range=[22:33];
   R_PA =[R_PA;R_TL(tile_range)'];
   Z_PA =[Z_PA;Z_TL(tile_range)'];
  DR_PA =[DR_PA;DR_TL(tile_range)'];
  DZ_PA =[DZ_PA;DZ_TL(tile_range)'];
  RES_PA=[RES_PA;RES_TL(tile_range)'];

 ff='00';
 for k=1:length(R_TL(tile_range));
  if(k>=10&k<99);ff='0';end;  if(k>=100);ff='';end;
  eval(['LA_PA=[LA_PA;''TL',ff,num2str(k,3),'F001''];']);
 end;




% DIAGNOSTICS
 R_BP=bpol_r;
 Z_BP=bpol_z;
 TH_BP=bpol_theta;  
 R_FL=f_r(1:38);
 Z_FL=f_z(1:38);


% for the bpol probes  
 LA_BP=[];
 for k=1:length(R_BP);
  ffk='00';
   if(k>=10&k<99);ffk='0';end;  if(k>=100);ffk='';end;
   eval(['LA_BP=[LA_BP;''BP',ffk,num2str(k,3),'F',ffk,num2str(k,3),'''];']);
  end;
% for the flux loops 
 LA_FL=[];
 for k=1:length(R_FL);
  ffk='00';
   if(k>=10&k<99);ffk='0';end;  if(k>=100);ffk='';end;
   eval(['LA_FL=[LA_FL;''FL',ffk,num2str(k,3),'F',ffk,num2str(k,3),'''];']);
  end;





clear I_T Ic Ifils Iphi 
clear a* b* c* d* e* f* g* h* i* j* k* l* m*
clear n* o* p* q* r* s* t* u* v* w* x* y* z*

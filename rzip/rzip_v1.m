%
%    RZIP - a linear tokamak plasma equilibrium response model
%    
%    1999 J Wainwright
%    
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
%
%    In addition, if you use RZIP code or a derivative, we ask that you reference
%    the appropriate papers both in any resulting code or publications. 
%
% A routine to create a filament Tokamak model. Cmatrix is via flux and fields.  Based on the rzip concept and methods as outlined in 
% the NF paper to be and originally coded by Idrinil
%
% Original version Dec 98
% updated to function call rzip and basic input/output list
% updated to expanded output/input and labelling convention changed
% updated to speed extensions, workspace checking etc
% updated to exact coupling from probes to plasma grid
% updated to be compatible with version2 standards. Full plasma grid, subroutines
%  
% John Wainwright Sep'99
%
% Inputs
% ======
% shot,         a reference to the equilibrium definition (will define I_AC, JG, RG, ZG, AG and BETA_P)
% time,         if reading real machine data the time at which to acquire the above
% used_coils,   select which coils are relevant to the experiment.
% n_eig,        how many eigenvalues to use.
% plas_ext,     [Rp, Rp', betap_dot, li_dot] externally defined plasma variables: plasma resistance, d(resistance)/dr
%               betapdot and lidot; the latter may or may not exist to define the disturbances.
% Machine_name, a flag for the machine configuration, must include the data outlined below
%               Machine_name is a reference to a local file and must be entered as characters.
%
%  NB The code produces 2 files depending on use. For a given machine it produces a file
%  ./Machine_name_plasma_less_v1.mat which has the purely structural details within it.
%  ./Machine_name_plasma_full_v1.mat which has the coupling from the structure to the proposed plasma grid.
%  Also due to the nature of the workspace checking, it is advised to 'clear all' before
%  creating a new model. Use lower case names to ensure VMS compatibility.
%
% active coils
% -------------
% R_AC,		radial position
% Z_AC.		vertical position
% RES_AC.	resistance
% DR_AC,	width
% DZ_AC,	height
% NT_AC,	number of turns in a coil
% I_AC,		equilibrium currents
% LA_AC,	label to distinguish coils: convention ???##F001
%                                                        ^  ^  ^
%                                                        |  |  L_ filament 1
%                                                        |  L____ Filament 
%                                                        L_______ name (5 chars 3/2 eg PF101)
% passive structure
% -----------------
% R_PA,		radial position
% Z_PA.		vertical position
% RES_PA.	resistance
% DR_PA,	width
% DZ_PA,	height
% LA_PA,	label to distinguish parts of structure: convention ??###F001
%                                                                     ^  ^  ^
%                                                                     |  |  L__ filament 1
%                                                                     |  L_____ Filament 
%                                                                     L________ name (5 chars 2/3 eg VV001)
% Diagnostics
% ----------- 
% R_BP,		radial position of Bpol probe
% Z_BP.		vertical position of Bpol probe
% TH_BP,        angle of Bpol probe
% R_FL,		radial position of flux probe
% Z_FL.		vertical position of flux probe
%
% plasma grid (grid over which plasma will be defined)
% ----------------------------------------------------
% RG, 		radial coordinate
% ZG,		vertical coordinate
% JG,		plasma current at point (RG,ZG) [If=0 then get plasmaless model]
% AG,		plasma current area at point (RG,ZG)  
% BETA_P,       an external input for beta_p consistant with above
% LI,           an external input for li consistant with above
%
%
% Outputs 
% =======
% A, B, C, D       for C as defined in the appropriate section.
%
% curlyM, curlyR   fundamental model matrix steps.
% Fs, Fd           stabilising and destabilising forces respectively
% decaybar, bz_bar, Current weighted fields terms.
% dbrdzbar, alpha
% IFLAG             A status flag for correct execution (0=fine)
%
% Code definition
% ===============
% n_eig_max        is a variable which defines the largest number of
%                  eigenmodes to be considered, memory useful
%

 function [A, B, C, D, curlyM, curlyR, Fs, Fd, Decaybar, bz_bar, dbrdzbar, alpha, IFLAG] = ...
                       rzip_v1(shot, time, used_coils, n_eig, plas_ext, Machine_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global variables associated with plasma grid
 global ZG_old RG_old  
 global dbzdz dbrdr dbzdr br bz bpg flg bprg flrg bpzg flzg mg             

% global variables associated with plasma less model
 global Le Re Er r z nr 
 global nac naf npf afri NT_AC
 global flcs bpcs

% define an upper limit to the number of eigenmodes to consider
  n_eig_max=40;
% flag for screen output on(1)/off(0)   
  screen_on=1;
% flag for whether interpolation or recalculation occurs yes(1)/no(0)
  do_interp=1;
% flag for whether code has executed properly (0 is true)
  IFLAG=0;

% Initialise situation dependent flags
% flag for decision on whether to make the plasma_less model (default no)
  make_plasma_less=0;
% flag for decision on whether to make the machine-plasmagrid links model (default no)
  make_plasma_full=0;
% flag for decision on whether to use workspace data, if it exists (default no)
  no_load=0;
% flag for decision on adding disturbance matrices (default no)
  make_disturb =0;
% flag for whether a plasma is required (default no)
  plasma_less = 0; 

% Create two file name variables to the output files.
  plasma_less_file  = [Machine_name,'_plasma_less_v1.mat'];
  plasma_full_file  = [Machine_name,'_plasma_full_v1.mat'];

% if not a local file check if is in the path structure somewhere.
  remote_plasma_less_file=which(plasma_less_file);
  remote_plasma_full_file=which(plasma_full_file);
  if(length(remote_plasma_less_file)); plasma_less_file=remote_plasma_less_file;  end;
  if(length(remote_plasma_full_file)); plasma_full_file=remote_plasma_full_file;  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE EXECUTION REQUIREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Is there a plasma (shot!=0)?
%=============================
 if(shot==0); 

% define the plasma_less flag
  plasma_less=1; 
  if(screen_on); disp('Assuming plasmaless'); end; 

 else; 

% get the plasma equilibrium details (I_AC, JG, ZG, RG, BETA_P, LI) 
% -----------------------------------------------------------------------
  [JG, AG, RG, ZG, BETA_P, LI, I_AC] = feval([Machine_name,'_plasma'],shot,time);
  ohm_p=plas_ext(1);   ohm_prime=plas_ext(2);   betap_dot=plas_ext(3);   li_dot=plas_ext(4);
  if(li_dot~=0 | betap_dot~=0) ; make_disturb=1 ;end;
% if we're using tcv then we should use all coils possible
  if(Machine_name(1:3)=='tcv');
   if(length(I_AC)<length(used_coils));used_coils=1:length(I_AC);end;
   if(max(used_coils)>19); used_coils=1:19; end;
   if(isnan(BETA_P)) ; IFLAG=1; break; end;
   if(isinf(BETA_P)) ; IFLAG=1; break; end;
  end;

 end;     % the plasma_less check loop

% Are the variables in the workspace (length(Le)>0) (NB don't use exist as initial global have size)
% will not reload unless i) plasma grids are different or ii) n_eig+nuc has changed (simplistic!) 
%===================================================================================================
 if(length(Le)>0);

% define the loading flag
  no_load=1;

% i) if the plasma grids are different interpolate or recalculate
% ---------------------------------------------------------------
% only check plasma grid if there is a plasma
%--------------------------------------------
  if(~plasma_less)   
   if(length(RG_old)~=length(RG) | length(ZG_old)~=length(ZG));
    make_plasma_full=1; 
   else;
    if(any(RG_old ~= RG) | any(ZG_old ~= ZG));
    make_plasma_full=1; 
    end;
   end;
  end;

% ii) if the number of coils or eigenvalues has changed, must reload (true for plasma and plasma_less cases)
% ----------------------------------------------------------------------------------------------------------
  if(length(Le)~=length(used_coils)+n_eig)
   if(screen_on);disp('The number of coils or eigenmodes has changed, must reload data');end;
   no_load=0;
  end

 else;  % variables are not in the workspace but maybe we have them on file

% check whether the plasmaless model has been already created  
%------------------------------------------------------------ 
  file_test_less = fopen(plasma_less_file,'r');
  if(file_test_less==-1); make_plasma_less=1; else; fclose(file_test_less); end; 
  

% check whether the coil-grid data model has already been created 
%---------------------------------------------------------------- 
  if(~plasma_less); 

   file_test_full=fopen(plasma_full_file,'r');;
   if(file_test_full==-1);
    make_plasma_full=1; 
   else;
% file probably exists but however this could be a new grid so check if we should recalculate
%--------------------------------------------------------------------------------------------
   fclose(file_test_full)             
   feval('load',plasma_full_file);
%  i) if the plasma grids are different  recalculate
%  ------------------------------------------------- 
    if(length(RG_old)~=length(RG) | length(ZG_old)~=length(ZG));
     make_plasma_full=1; 
    else;
     if(any(RG_old ~= RG) & any(ZG_old ~= ZG));
     make_plasma_full=1; 
     end;
    end;
   end;     
  end;   % end of plasma_less check                  
 end;    % end the check for variables in the workspace  

 if(~plasma_less);   NG = length(RG); end;

% start clock
%------------
 tic 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there is no plasma less model we will have to create it %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if(make_plasma_less) 

  t_start_less=toc;
  f_start_less=flops;

% tell the user what's happening
%-------------------------------- 
  if(screen_on);disp('Found no saved plasmaless data; creating from scratch');end;

% get the relevant machine details
%--------------------------------- 
  eval(['[R_AC,Z_AC,RES_AC,DR_AC,DZ_AC,LA_AC,NT_AC,R_PA,Z_PA,RES_PA,DR_PA,DZ_PA,LA_PA,R_BP,Z_BP,TH_BP,LA_BP,R_FL,Z_FL,LA_FL]=',Machine_name,';']);

% To ensure consistancy work only with column inputs so make sure all arrays are columns;
%-----------------------------------------------------------------------------------------
   R_AC=check_is_column(R_AC);   Z_AC=check_is_column(Z_AC); RES_AC=check_is_column(RES_AC);
  DR_AC=check_is_column(DR_AC); DZ_AC=check_is_column(DZ_AC); LA_AC=check_is_column(LA_AC);  NT_AC=check_is_column(NT_AC);

   R_PA=check_is_column(R_PA);   Z_PA=check_is_column(Z_PA); RES_PA=check_is_column(RES_PA);
  DR_PA=check_is_column(DR_PA); DZ_PA=check_is_column(DZ_PA); LA_PA=check_is_column(LA_PA);

   R_BP=check_is_column(R_BP);   Z_BP=check_is_column(Z_BP);  TH_BP=check_is_column(TH_BP);  LA_BP=check_is_column(LA_BP);
   R_FL=check_is_column(R_FL);   Z_FL=check_is_column(Z_FL);  LA_FL=check_is_column(LA_FL);

% First of all we need to make decisions about how many actual coils there are
% as opposed to circuits. This is information stored in LA (for label);
% initially we must work with filaments and combine to the coils.
%------------------------------------------------------------------------------

% active coils (nac = number of active coils, naf = number of active filaments)
%==============------------------------------------------------------------------
% seek out the different coils by name
 [nac,naf,afri,coil_names]=unique_name(LA_AC(:,1:3));

% passive coils (nps = number of passive 'coils', npf = number of passive filaments)
%===============---------------------------------------------------------------------
% seek out the different structures by name
 [nps,npf,pfri,pass_names]=unique_name(LA_PA(:,1:2));

 if(screen_on);
  disp(['Potentially ',num2str(nac),' Active coils and ',num2str(nps),' Passive structure types']);
  fill=num2str(length(used_coils));if(length(used_coils)==nac);fill='all of the';end;
  disp(['will use ',fill,' Active coils and ',num2str(n_eig),' eigenmodes']);
  disp(['There are ',num2str(length(R_BP)),' Bpol probes and ',num2str(length(R_FL)),' flux loops']);
 end;

% however will treat the separate passive structures as one structure (plotting use).
 nps=1;pfri(1,1)=1;prfi(1,2)=npf;

% so we can begin the inductance matrix calculation
%===================================================

% from the active filaments to the passive filaments
  for k=1:naf
   m_af_pf(k,1:npf)=mutual_ind(R_AC(k),R_PA,Z_PA-Z_AC(k))';
  end;

% from the active filaments to the active filaments
  for k=1:naf
   m_af_af(k,1:naf)=mutual_ind(R_AC(k),R_AC,Z_AC-Z_AC(k))';
  end;

% from the passive filaments to the passive filaments
  for k=1:npf
   m_pf_pf(k,1:npf)=mutual_ind(R_PA(k),R_PA,Z_PA-Z_PA(k))';
  end;


% obviously there will be some difficulties we need to
% redefine the self inductances and to make sure that if
% coils are too close we are accurate enough. (should be ok here)
% Also remember mutual_ind removes NaN and Inf and replaces with 0

% m_af_af=m_af_af+diag(self_circ(R_AC,sqrt((DR_AC.*DZ_AC)/pi)));
 m_af_af=m_af_af+diag(self_sqre(R_AC,DR_AC,DZ_AC));
 m_pf_pf=m_pf_pf+diag(self_sqre(R_PA,DR_PA,DZ_PA));

% need to combine filaments to coils and include the number of turns  
%-------------------------------------------------------------------

% for the active/passive connections
  M_a_p= coil_from_filaments(m_af_pf.*[NT_AC*ones(1,npf)],afri,nac,naf,0);

% for the active connections
% compact rows
 temp = coil_from_filaments(m_af_af.*[NT_AC*ones(1,naf)],afri,nac,naf,0)';
% then columns  
 L_a_a= coil_from_filaments(temp.*[NT_AC*ones(1,nac)],afri,nac,naf,0)';

% Memory saving
 clear m_af_af m_af_pf

 if(screen_on);disp(['Created the plasmaless L and R after ',num2str(toc-t_start_less),'s']);end;

% assume there is no compacting for the VV, if there is it would
% be achieved using a NT_PA variable

% construct the eigenmodes (see notes) 
%=====================================

% define eigenfunctions (vectors Ee and values Ve)
 inv_R  =inv(diag(RES_PA));
 inv_R_L=inv_R*m_pf_pf;
 [Ee,Ve]=eig(inv_R_L);

% sort (V_o) into ascending importance of values (vectors Es and values Vs)
 [i,V_o]=esort(diag(Ve));
 Es=Ee(:,V_o); 
 temp = (diag(Ve)); temp= temp(V_o) ;  Vs=diag(temp); clear temp;

% renormalise to net passive resistance
% eigencurrents will be eigenmodes but scaled by norm_mat
% normalising matrix is diagonal and equal to transpose.
  norm_res=1/sum(1./RES_PA);
  norm_mat=sqrt(norm_res*inv(Es)*inv_R*inv(Es'));
  norm_mat=diag(diag(norm_mat));

% rescale the eigenvectors appropriately
  Es=Es*norm_mat;
  Vs=Vs*norm_res;

% set maximum number of emodes
  if(n_eig_max>npf)
   n_eig_max=npf; if(screen_on);disp('More eigenmodes than points');end;
  end;  

% define the reduced set of eigenfunctions (vectors Er and values Vr)
  Er=Es(:,1:n_eig_max); 
  Vr=Vs(1:n_eig_max,1:n_eig_max); 

% define the equivalent inductance matrix
%----------------------------------------
  Le=[L_a_a         M_a_p*Er;
     [M_a_p*Er]'        Vr ];

% and the equivalent Resistance matrix 
%--------------------------------------
  Re=[diag(coil_from_filaments(RES_AC,afri,nac,naf,0)) zeros(nac,n_eig_max);zeros(n_eig_max,nac) norm_res*eye(size(Vr))]; 

% Memory saving
  clear m_pf_pf M_a_p inv_R inv_R_L Ee Es Ve V_o Vr Vs 

 if(screen_on);disp(['Eigenmode adjusted L and R after ',num2str(toc-t_start_less),'s']);end;

% define vectors of filament positions (coils and structure)
  r=[R_AC;R_PA];
  z=[Z_AC;Z_PA];
  nr=length(r);

% C Matrix information takes some construction (we're considering b-pol and flux loops)
%======================================================================================
 
  nbp = length(R_BP); nfl = length(R_FL); 

% Consider coupling of bpol probes (bp) to coils and structure (cs)
%------------------------------------------------------------------
  Bprcs=zeros(nbp,nr);Bpzcs=zeros(nbp,nr);
  for k=1:nbp;
   Bprcs(k,:)  =   bfield_br(r,R_BP(k),Z_BP(k)-z,1)'; 
   Bpzcs(k,:)  =   bfield_bz(r,R_BP(k),Z_BP(k)-z,1)'; 
  end;

% convert to eigenmodes
  Bprcs = [Bprcs(:,1:naf) Bprcs(:,naf+1:naf+npf)*Er];
  Bpzcs = [Bpzcs(:,1:naf) Bpzcs(:,naf+1:naf+npf)*Er];
 
% but have to create the net effects of the links to the coils

  Bprcs  = Bprcs'.*[[NT_AC;ones(n_eig_max,1)]*ones(1,nbp)];
  Bpzcs  = Bpzcs'.*[[NT_AC;ones(n_eig_max,1)]*ones(1,nbp)]; 
  bprcs   = coil_from_filaments(Bprcs,afri,nac,naf,n_eig_max)';
  bpzcs   = coil_from_filaments(Bpzcs,afri,nac,naf,n_eig_max)';

% Create total field(bpt) and angle(bpa)
  bptcs   = sqrt(bprcs.*bprcs+bpzcs.*bpzcs); 
  bpacs   = atan2(bpzcs,bprcs); 

% compare with angle of coil
  net_angle=bpacs-[TH_BP*ones(1,nac+n_eig_max)];

% so the final effect from coils/structure to bpols is
  bpcs=cos(net_angle).*bptcs;


% and now the flux probes (fl) (which is simply via mutual inductance F=LI)
%---------------------------------------------------------------------------
% first to the coils and structure (cs)
  Flcs=zeros(nfl,nr);
  for k=1:nfl;
   Flcs(k,:)  =     mutual_ind(r,R_FL(k),Z_FL(k)-z)'; 
  end;
% convert to eigenmodes
  Flcs = [Flcs(:,1:naf) Flcs(:,naf+1:naf+npf)*Er];
  
% but have to create the net effects of the links to the coils
  Flcs   = Flcs'.*[[NT_AC;ones(n_eig_max,1)]*ones(1,nfl)];
  flcs   = coil_from_filaments(Flcs,afri,nac,naf,n_eig_max)';

% Save the necessary outputs for later
  old_clock=clock;
  feval('save', plasma_less_file, ...
               'n_eig_max','npf','nac','naf','afri','NT_AC','r','z','nr','bpcs','flcs','Le','Re','Er','old_clock');
  if(screen_on);
   disp(['Creating the plasmaless model required ',num2str(toc-t_start_less),'s and ',num2str(flops-f_start_less,8),' flops']);
  end;

  clear Flcs bpacs bptcs bprcs bpzcs net_angle Bprcs Bpzcs norm_mat

 else;   % read in the old plasma_less data
         %---------------------------------

% tell the user what's happening
%-------------------------------- 
  if(no_load)
   if(screen_on);disp('Using workspace plasma less data');end;
  else
   feval('load',plasma_less_file);
   if(screen_on);
    cre_date=[num2str(old_clock(3)),'/',num2str(old_clock(2)),'/',num2str(old_clock(1))];
    if(old_clock(4)<10);fil1='0';else;fil1='';end;
    if(old_clock(5)<10);fil2=':0';else;fil2=':';end;
    cre_time=[fil1,num2str(old_clock(4)),fil2,num2str(old_clock(5))];
    disp(['Using saved plasmaless data from ',cre_time,' on the ',cre_date])
   end;
% check the current number of eigenmodes is less than the max
  if(n_eig>n_eig_max)
   n_eig=n_eig_max;
   if(screen_on);disp(['n_eig larger than this saved data allows, reset to ',num2str(n_eig_max)]);end;
   end;
  end;

 end;   % make_plasma_less loop
        %----------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the appropriate coils and eigenmodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the number of used coils (nuc) and define the matrix used_currents
% -------------------------------------------------------------------------
 if(no_load)
  nuc=length(Le)-n_eig;
  used_currs=[1:(nuc+n_eig)];
 else
  used_currs=[used_coils,nac+1:nac+n_eig];
  nuc=length(used_coils);
  Le=Le(used_currs,used_currs);
  Re=Re(used_currs,used_currs);
  bpcs=bpcs(:,used_currs);
  flcs=flcs(:,used_currs);
 end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there are no links from the machine to the grid the will have to create them %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if(make_plasma_full) 

  t_start_full=toc;
  f_start_full=flops;

% tell the user what's happening
%-------------------------------- 
 if(screen_on);disp('Found no saved plasmafull data; creating from scratch');end;

% however at this stage if we haven't performed the plasma_less
% model then we need to reobtain the machine data
%------------------------------------------------
 if(~make_plasma_less) ;
  eval(['[R_AC,Z_AC,RES_AC,DR_AC,DZ_AC,LA_AC,NT_AC,R_PA,Z_PA,RES_PA,DR_PA,DZ_PA,LA_PA,R_BP,Z_BP,TH_BP,LA_BP,R_FL,Z_FL,LA_FL]=',Machine_name,';']);
  nbp = length(R_BP); nfl = length(R_FL); 
  r=[R_AC;R_PA];  z=[Z_AC;Z_PA];  nr=length(r);
 end;

% we need to define various interactions of the structure to the plasma grid. Assume the number of grid points is the larger
% and therefore optimise the calculations to this. We need Bz, Br, Mut(a/p,g) dBzdr (the bulk of the code works in this bit) 

  TMP=zeros(nr,NG);
  for k=1:nr;
   TMP(k,:) = bfield_bz(r(k),RG,ZG-z(k),1)'; 
  end;

% with an eigenmode analysis we need to adjust these fields to be scaled by the eigenvectors.
% Ie we convert [X_af X_ps]*[I_af I_ps] to [X]*[I_af I_e]  X=[X_af X_ps*E] 
% (bear in mind the transposition due to definition of eg bz to grid)
   TMP = [TMP(1:naf,:);[TMP(naf+1:naf+npf,:)'*Er]'];

% but have to create the net effects of the links to the coils
% define a Number of turns array
   TMP=[[NT_AC;ones(n_eig_max,1)]*ones(1,NG)].*TMP;
   bz = coil_from_filaments(TMP,afri,nac,naf,n_eig_max);
   if(screen_on);disp('Finished with bz');end;

   TMP=zeros(nr,NG);
   for k=1:nr;
    TMP(k,:) = bfield_br(r(k),RG,ZG-z(k),1)'; 
   end;
   TMP = [TMP(1:naf,:);[TMP(naf+1:naf+npf,:)'*Er]'];
   TMP=[[NT_AC;ones(n_eig_max,1)]*ones(1,NG)].*TMP;
   br    = coil_from_filaments(TMP,afri,nac,naf,n_eig_max);
   if(screen_on);disp('Finished with br');end;

   TMP=zeros(nr,NG);
   for k=1:nr;
    TMP(k,:) = mutual_ind(r(k),RG,ZG-z(k))'; 
   end;
   TMP = [TMP(1:naf,:);[TMP(naf+1:naf+npf,:)'*Er]'];
   TMP=[[NT_AC;ones(n_eig_max,1)]*ones(1,NG)].*TMP;
   mg    = coil_from_filaments(TMP,afri,nac,naf,n_eig_max);
   if(screen_on);disp('Finished with mg');end;

   TMP=zeros(nr,NG);
   for k=1:nr;
    TMP(k,:) = bfield_dbzdr(r(k),RG,ZG-z(k),1)'; 
   end;
   TMP = [TMP(1:naf,:);[TMP(naf+1:naf+npf,:)'*Er]'];
   TMP=[[NT_AC;ones(n_eig_max,1)]*ones(1,NG)].*TMP;
   dbzdr = coil_from_filaments(TMP,afri,nac,naf,n_eig_max);
   clear TMP Er
   if(screen_on);
    disp('Finished with dbzdr');
    disp(['Creating the main matrices took ',num2str(toc-t_start_full),' seconds']);
   end;

% plasma grid mapping to the diagnostics. Using exact fields and using new routines 
%==================================================================================

% first get coupling of bpol probes to plasma grid (g)
%-----------------------------------------------------
  bprg=zeros(nbp,NG);bpzg=zeros(nbp,NG);
  bpdrdrg=zeros(nbp,NG);bpdzdrg=zeros(nbp,NG);bpdzdzg=zeros(nbp,NG);
  for k=1:nbp;
   bprg(k,:)     =    bfield_br(RG,R_BP(k),Z_BP(k)-ZG,1)'; 
   bpzg(k,:)     =    bfield_bz(RG,R_BP(k),Z_BP(k)-ZG,1)'; 
   bpdzdrg(k,:)  =    bfield_dbzda(RG,R_BP(k),Z_BP(k)-ZG,1)'; 
   bpdzdzg(k,:)  =    bfield_dbzdc(RG,R_BP(k),Z_BP(k)-ZG,1)'; 
   bpdrdrg(k,:)  =    bfield_dbrda(RG,R_BP(k),Z_BP(k)-ZG,1)'; 
   bpdrdzg(k,:)  =    bfield_dbrdc(RG,R_BP(k),Z_BP(k)-ZG,1)'; 
  end;

% Create total field(bpt) and angle(bpa)
  bptg   = sqrt(bprg.*bprg+bpzg.*bpzg); 
  bpag   = atan2(bpzg,bprg); 

% create the linear terms (see the algebra)
  dthetadr = (bprg.*bpdzdrg-bpzg.*bpdrdrg)./bptg;
  dthetadz = (bprg.*bpdzdzg-bpzg.*bpdrdzg)./bptg;
  dmagdr   = (bprg.*bpdrdrg+bpzg.*bpdzdrg)./bptg;
  dmagdz   = (bprg.*bpdrdzg+bpzg.*bpdzdzg)./bptg;

% compare with angle of coil
  net_angle=[TH_BP*ones(1,NG)]-bpag;

% so the final effect from grid to bpols is
  bpg=cos(net_angle).*bptg;
  bprg=sin(net_angle).*dthetadr+cos(net_angle).*dmagdr;
  bpzg=sin(net_angle).*dthetadz+cos(net_angle).*dmagdz;

% and now the flux probes (fl) (which is simply via mutual inductance F=LI)
%--------------------------------------------------------------------------
  flg=zeros(nfl,NG);  VPr=zeros(nfl,NG); VPz=zeros(nfl,NG);
  for k=1:nfl;
   flg(k,:)  =     mutual_ind(RG,R_FL(k),Z_FL(k)-ZG)'; 
   VPr(k,:)  =     vecpot_dada(RG,R_FL(k),Z_FL(k)-ZG,1)'; 
   VPz(k,:)  =     vecpot_dadc(RG,R_FL(k),Z_FL(k)-ZG,1)'; 
  end;

  flrg  =  2*pi*[R_FL*ones(1,NG)].*VPr;
  flzg  =  2*pi*[R_FL*ones(1,NG)].*VPz;

  RG_old=RG; ZG_old=ZG; old_clock=clock;
  feval('save',plasma_full_file, ...
              'flg','flrg','flzg','bpg','bprg','bpzg','bz','br','dbzdr','mg','old_clock','RG_old','ZG_old');
  if(screen_on);
   disp(['Creating the links to the plasma grid required ',num2str(toc-t_start_full),'s and ',num2str(flops-f_start_full,8),' flops']);
  end;

 else;   % read in the old plasma_full data

% tell the user what's happening
%-------------------------------- 
% if no_load is switched on it may be just due to the plasmaless considerations; therefore check whether bz is really in the workspace.
  if(no_load & length(bz)>0)
   if(screen_on);disp('Using workspace plasma full data');end;
  else
   if(~plasma_less)
    feval('load',plasma_full_file);
    if(screen_on);
     cre_date=[num2str(old_clock(3)),'/',num2str(old_clock(2)),'/',num2str(old_clock(1))];
     if(old_clock(4)<10);fil1='0';else;fil1='';end;
     if(old_clock(5)<10);fil2=':0';else;fil2=':';end;
     cre_time=[fil1,num2str(old_clock(4)),fil2,num2str(old_clock(5))];
     disp(['Using saved plasmafull data from ',cre_time,' on the ',cre_date])
    end;
   end;
  end;

 end;   % make_plasma_full


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the appropriate coils and eigenmodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if(~plasma_less);
  bz=bz(used_currs,:);
  br=br(used_currs,:);
  dbzdr=dbzdr(used_currs,:);
  mg=mg(used_currs,:);
 end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make rzip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% And then make a proper plasmaless or rzip model
%===============================================

 if(plasma_less)
  curlyM = Le ;
  curlyR = Re ;
  Cbpol= bpcs;
  Cpsi= flcs;
  Cipol=[eye(nuc) zeros(nuc,n_eig)];

% set the plasma output terms to nonexistent

  Fs=[];Fd=[];Decaybar=[];bz_bar=[];dbrdzbar=[];alpha=[];

 else; %%%%%%%%%%%%% else of plasma less check

% define mu0
  mu0=4e-7*pi;

% define the equilibrium circuit currents
 i_full=[I_AC(used_coils);zeros(n_eig,1)];

% plasma current ~ int JG dA
 DA=mean(AG);
 IP=DA*sum(JG);
 J_SUM=sum(JG); 

% need to have read in Beta_p  externally
 betap=BETA_P;
 li=LI;

% plasma mutual inductance gradients
%-----------------------------------
 mrprime = 2*pi*(bz*[RG.*JG])/J_SUM;
 mzprime =-2*pi*(br*[RG.*JG])/J_SUM;

% plasma field gradients
%----------------------- 
 dbzdr_plas = i_full'*dbzdr;
 dbzdrbar   = (dbzdr_plas*JG)/J_SUM;
 dbrdzbar   =  dbzdrbar;
 rdbrdzbar  = (dbzdr_plas*[JG.*RG])/J_SUM;

% Useful terms for tokamak people
 alpha =  2*pi*rdbrdzbar/IP;
 bz_x = bz'*i_full; 
 bz_bar= bz_x'*JG/J_SUM;
 rbz_bar= bz_x'*(JG.*RG)/J_SUM;
 Decaybar = -dbzdr_plas*(JG.*RG./bz_x)/J_SUM;
 Fd = -2*pi*IP*rdbrdzbar;
 ind_pa=nuc+1:nuc+n_eig;
 Fs = mzprime(ind_pa)'*(mzprime(ind_pa)./diag(Le(ind_pa,ind_pa)))*IP*IP; 

% plasma properties
%------------------  
 R0    =  sum(RG'*JG)/J_SUM; 
 gamma0 = -4*pi*R0*bz_bar/mu0/IP; fo = gamma0-betap-1/2;
 Lp0 = mu0*R0*fo;
 mut_bar = mg*JG/J_SUM;


% display useful coefficients
%----------------------------
 if(screen_on);disp(['bz_bar,alpha = ', num2str(bz_bar),'  ',num2str(alpha)]);end;


% create the model 'inductance' and 'resistance' matrices. Is a reduced model as ZIp has been elimated.
%======================================================================================================
 M11 = Le+mzprime*mzprime'/alpha;  M12 = mrprime*IP;                                M13 = mut_bar; 
 M21 = M12';                       M22 = (mu0*IP/2/R0+2*pi*(bz_bar+rdbrdzbar))*IP;  M23 = mu0*IP*gamma0+2*pi*rbz_bar;
 M31 = mut_bar';                   M32 = mu0*(1+fo)*IP+2*pi*rbz_bar;                M33 = Lp0;

 curlyM = [ M11  M12  M13 ; 
            M21  M22  M23 ; 
            M31  M32  M33 ];

 curlyR = [ Re  zeros(length(Re),2);
            zeros(1,length(curlyM));
            zeros(1,length(curlyM)-2) ohm_prime*IP ohm_p];   


% rigid plasma movement contribution to bpol probes 
  bp_bar = bpg*JG/sum(JG);
  dbpdr  = bprg*JG/sum(JG)*IP;
  dbpdz  = bpzg*JG/sum(JG);
% rigid plasma movement contribution to flux loops 
  fl_bar = flg*JG/sum(JG);
  dfldr  = flrg*JG/sum(JG)*IP;
  dfldz  = flzg*JG/sum(JG);

  Cr   = [zeros(1,nuc+n_eig) 1 0];
  Cip  = [zeros(1,nuc+n_eig) 0 1];
  Czip = [mzprime'/alpha 0 0];
  Cipol= [eye(nuc) zeros(nuc,n_eig+2)];
  Cbpol= [bpcs dbpdr bp_bar] + dbpdz*Czip;
  Cpsi = [flcs dfldr fl_bar] + dfldz*Czip;

 end; %%%%%%%%%%%%% end of plasma less check


% Make the A B C and D matrices
% -----------------------------
 minv = inv(curlyM);
 extra_zero=0;if(~plasma_less);extra_zero=2;end;
 volt = [diag(ones(nuc,1)); zeros(n_eig+extra_zero,nuc)]; 

 A =-minv*curlyR;
 B = minv*volt;
 C  = [Cbpol; Cpsi; Cipol];
 D  = zeros(size(C,1),size(B,2));

% If we have disturbances then we need to rearrange things
 if(make_disturb)
  e_dist=-mu0*IP*IP/2;
  E=minv*[zeros(nuc+n_eig,2); e_dist e_dist/2; 0 0];
% Rearange states so that Z=x-Ew ; and V=[v_coils beta_dot li_dot];
  B=[B A*E];
  D=[D C*E];
 end;

 if(screen_on);
  eA=esort(eig(A));disp('First 10 eigenmodes are');
  for k=1:10;disp(num2str(eA(k)));end;
 end;

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES required by RZIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       bfield_br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the radial b field due to a circular filament.
% a is radius of filament, b is radius of position,  c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function br = bfield_br(a,b,c,cur)

% set up constants
   mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;

% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;

  f=com*1./sqrt((a+b).^2+c.^2);
  f1=-K+(a.*a+b.*b+c.*c)./((a-b).^2.+c.^2.).*E;

% define the value of br (Smythe pg 290)
  br=2.e-7*cur.*c.*f./b.*f1;

  clear E K f f1 k2 mu0;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          bfield_bz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the vertical b field due to a circular filament.
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function bz = bfield_bz(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;

% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
 [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;

  f=com*1./sqrt((a+b).^2.+c.^2.);
  f2=K+(a.*a-b.*b-c.*c)./((a-b).^2.+c.^2.).*E;

% define the value of bz (Smythe pg 290)

  bz=2.e-7.*cur.*f.*f2;

  clear f f2 mu0 k2 E K;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             dbrda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the rate of change of the radial magnetic field of a current loop at a fixed point due to 
% radial movement of the source.
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbrda = bfield_dbrda(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;


% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
  [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;
  dEdk= (E./k-K./k);
  dKdk= (E./((1-k.^2).*k)-K./k);

  r=b;
  z=c;
  ak=k;

     g2=(a-r).^2+z.^2;
  dg2da=2*(a-r);

     g3=(a+r).^2+z.^2;
  dg3da=2*(a+r);

     g4=a.^2+r.^2+z.^2;
  dg4da=2*a;

  dkda=(4*r.*g3-4*dg3da.*a.*r)./(2*k.*g3.^2);


  F3=-K+g4.*E./g2;
  F4=g3.^(-1/2);

  dF3da=-dKdk.*dkda+(dg4da.*E+dEdk.*dkda.*g4)./g2 - dg2da.*g4.*E./(g2.^2);
  dF4da=-dg3da./(2*g3.^(3/2));

  dbrda= 2e-7*cur.*z./r.*(dF3da.*F4+dF4da.*F3);

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             dbrdc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the rate of change of the radial magnetic field of a current loop at a fixed point due to 
% vertical movement of the source.
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbrdc = bfield_dbrdc(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;


% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
  [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;
  dEdk= (E./k-K./k);
  dKdk= (E./((1-k.^2).*k)-K./k);

  r=b;
  z=c;
  ak=k;

     g2=(a-r).^2+z.^2;
  dg2dz=2*z;

     g3=(a+r).^2+z.^2;
  dg3dz=2*z;

     g4=a.^2+r.^2+z.^2;
  dg4dz=2*z;

  dkdz=-(2*a.*r.*dg3dz)./(k.*g3.^2);


  F3=-K+g4.*E./g2;
  F4=g3.^(-1/2);

  dF3dz=-dKdk.*dkdz+(dg4dz.*E+dEdk.*dkdz.*g4)./g2 - dg2dz.*g4.*E./(g2.^2);
  dF4dz=-dg3dz./(2*g3.^(3/2));

  dbrdc=-2e-7*cur./r.*(dF3dz.*F4.*z+dF4dz.*F3.*z+F3.*F4);

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     bfield_dbrdr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the radial gradient of the radial magnetic field due to a circular filament
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbrdr = bfield_dbrdr(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;


% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
 [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;

  r=b;
  z=c;
  ak=k;

  g1=a.*a+r.*r+z.*z;
  g1r=2.*r;

  g2=(a-r).*(a-r)+z.*z;
  g2r=2.*r-2.*a;

  g3=(a+r).*(a+r)+z.*z;
  g3r=2.*r+2.*a;

  akr=(-2.*a.*r.*(a+r)./(g3.*g3)+a./g3)./sqrt(a.*r./g3);
  Er=E./ak-K./ak;
  Kr=E./(ak.*(com*1.-ak.*ak))-K./ak;

  p1=-2.0e-7*z.*(-K+E.*g1./g2)./(r.*r.*sqrt(g3));
  p2=-1.0e-7*z.*(-K+E.*g1./g2).*g3r./(r.*g3.^(1.5));
  p3= 2.0e-7*z.*(E.*g1r./g2-E.*g1.*g2r./(g2.*g2)+(g1.*Er.*akr)./(g2)-(Kr.*akr))./(r.*sqrt(g3));

  dbrdr=(p1+p2+p3).*cur;

   clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          bfield_dbzda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the rate of change of the vertical magnetic field of a current loop at a fixed point due 
% to radial movement of the source.
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbzda = bfield_dbzda(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;


% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
  [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;
  dEdk= (E./k-K./k);
  dKdk= (E./((1-k.^2).*k)-K./k);

  r=b;
  z=c;
  ak=k;


     g1=a.^2-r.^2-z.^2;
  dg1da=2*a;

     g2=(a-r).^2+z.^2;
  dg2da=2*(a-r);

     g3=(a+r).^2+z.^2;
  dg3da=2*(a+r);

  dkda=(4*r.*g3-4*dg3da.*a.*r)./(2*k.*g3.^2); 

  F1=K+g1.*E./g2;
  F2=g3.^(-1/2);

  dF1da=dKdk.*dkda+(dg1da.*E+dEdk.*dkda.*g1)./g2 - dg2da.*g1.*E./(g2.^2);
  dF2da=-dg3da./(2*g3.^(3/2));

  dbzda= 2e-7*cur.*(dF1da.*F2+dF2da.*F1);

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        bfield_dbzdc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the rate of change of the vertical magnetic field of a current loop at a fixed point due 
% to radial movement of the source.
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbzdc = bfield_dbzdc(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;

% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
  [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;
  dEdk= (E./k-K./k);
  dKdk= (E./((1-k.^2).*k)-K./k);

  r=b;
  z=c;
  ak=k;

     g1=a.^2-r.^2-z.^2;
  dg1dz=-2*z;

     g2=(a-r).^2+z.^2;
  dg2dz=2*z;

     g3=(a+r).^2+z.^2;
  dg3dz=2*z;

  dkdz=-(2*a.*r.*dg3dz)./(k.*g3.^2);

  F1=K+g1.*E./g2;
  F2=g3.^(-1/2);
  dF1dz=dKdk.*dkdz+(dg1dz.*E+dEdk.*dkdz.*g1)./g2 - dg2dz.*g1.*E./(g2.^2);   
  dF2dz=-dg3dz./(2*g3.^(3/2));

  dbzdc=-2e-7*cur.*(dF1dz.*F2+dF2dz.*F1);

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    bfield_dbzdr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the radial gradient of the vertical magnetic field due to a circular filament.
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbzdr = bfield_dbzdr(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;


% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
 [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;

  r=b;
  z=c;
  ak=k;

   g1=a.*a-r.*r-z.*z;
  g1r=-2.*r;

   g2=(a-r).*(a-r)+z.*z;
  g2r=2.*r-2.*a;

   g3=(a+r).*(a+r)+z.*z;
  g3r=2.*r+2.*a;

  akr=(-2.*a.*r.*(a+r)./(g3.*g3)+a./g3)./sqrt(a.*r./g3);
  Er=E./ak-K./ak;
  Kr=E./(ak.*(com*1.-ak.*ak))-K./ak;

  p1=-1.0e-7*(K+E.*g1./g2).*g3r./g3.^(3./2.);
  p2= 2.0e-7*(E.*g1r./g2-E.*g1.*g2r./(g2.*g2)+(g1.*Er.*akr)./g2+Kr.*akr)./sqrt(g3);

  dbzdr=(p1+p2).*cur;

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        bfield_dbzdz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the vertical gradient of the vertical magnetic field due to a circular filament
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dbzdz = bfield_dbzdz(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;


% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
  [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;

  r=b;
  z=c;
  ak=k;

   g1=a.*a-r.*r-z.*z;
  g1z=-2.*z;

   g2=(a-r).*(a-r)+z.*z;
  g2z=2.*z;

   g3=(a+r).*(a+r)+z.*z;
  g3z=2.*z;

  akz=(-a.*r.*g3z./(g3.*g3))./sqrt(a.*r./g3);

  Er=E./ak-K./ak;
  Kr=E./(ak.*(com*1.-ak.*ak))-K./ak;

  p1=-1.0e-7*(K+E.*g1./g2).*g3z./g3.^(3./2.);
  p2= 2.0e-7*(E.*g1z./g2-E.*g1.*g2z./(g2.*g2)+(g1.*Er.*akz)./(g2)+(Kr.*akz))./sqrt(g3);

  dbzdz=(p1+p2).*cur;

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2

  return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      check_is_column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a function to ensure vectors are column oriented;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [A]=check_is_column(A);

% only correct vectors
 if(any(size(A)==1))
  if(size(A,1)<size(A,2));A=A';end;
 end;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    coil_from_filaments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% link filaments to coils in the appropriate manner
% old_data and new_data are in order: coils then VV. 
% nc is number of coils. nf1 is number of filaments for the coils
% nf2 is the number of filaments for the VV, link  is the coil linkage information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [new_data]= coil_from_filaments(old_data,link,nc,nf1,nf2);

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
%                                        plot_box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a box of dimensions dr, dz at r,z in colour col.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function plot_box(r,z,dr,dz,col);

  nr=length(r);
  for k=1:nr 
   x=[r(k)-dr(k) r(k)-dr(k) r(k)+dr(k) r(k)+dr(k) r(k)-dr(k)];
   y=[z(k)+dz(k) z(k)-dz(k) z(k)-dz(k) z(k)+dz(k) z(k)+dz(k)];
   fill(x,y,col);
  end;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         plot_box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a box of dimensions dr, dz at r,z in faceclour colf and border colour colb.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function plot_box2(r,z,dr,dz,colb,colf);

 nr=length(r);
 for k=1:nr 
  x=[r(k)-dr(k) r(k)-dr(k) r(k)+dr(k) r(k)+dr(k) r(k)-dr(k)];
  y=[z(k)+dz(k) z(k)-dz(k) z(k)-dz(k) z(k)+dz(k) z(k)+dz(k)];
  fill(x,y,colf);
  x=[r(k)-dr(k) r(k)-dr(k) r(k)+dr(k) r(k)+dr(k) r(k)-dr(k)];
  y=[z(k)+dz(k) z(k)-dz(k) z(k)-dz(k) z(k)+dz(k) z(k)+dz(k)];
  plot(x,y,colb);
 end;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       self_circ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the self inductance of a coil of circular cross section, radius r and major radius a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [sel_ind] = self_circ(a,r)

 if(size(a,1)  <size(a,2))  ; a=a'; end;
 if(size(r,1)  <size(r,2))  ; r=r'; end;

 mu0=4e-07*pi;

 sel_ind = mu0*a.*(log(8*a./r)-1.75);

 clear aphi mu0 k;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       unique_name 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a file which give a colum vector of text names, counts the different names(nn) and defines the
% range in the list of each name (index range) and total length of string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            vecpot_dada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define rate of change of the vector potential at a point due to radial movement of a circular
% source. 
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dAda = vecpot_dada(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;

% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
  [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;
  dEdk= (E./k-K./k);
  dKdk= (E./((1-k.^2).*k)-K./k);

  r=b;
  z=c;
  ak=k;

     g3=(a+r).^2+z.^2;
  dg3da=2*(a+r);

  dkda=(4*r.*g3-4*dg3da.*a.*r)./(2*k.*g3.^2);

  F3=(2-k2).*K-2*E;
  dF3da=(2-k2).*dKdk.*dkda - 2*k.*dkda.*K - 2*dEdk.*dkda;

  dAda= 2e-7*cur./sqrt(r).*( sqrt(a).*dF3da./k + F3./(2*sqrt(a).*k) - sqrt(a).*dkda.*F3./k2 );

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         vecpot_dadc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define rate of change of the vector potential at a point due to vertical movement of a circular
% source. 
% a is radius of filament, b is radius of position, c is separation and cur is filament current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dAdc = vecpot_dadc(a,b,c,cur)

% set up constants
  mu0=4e-07*pi;

% ensure matrices are the same orientation
  if(size(a,1)  <size(a,2))  ; a=a'; end;
  if(size(b,1)  <size(b,2))  ; b=b'; end;
  if(size(c,1)  <size(c,2))  ; c=c'; end;
  if(size(cur,1)<size(cur,2)); cur=cur'; end;

% set up a dummy mapping array
  com=ones(size(a)); if(length(cur)==1) ; cur=cur*com; end;

% construct the elliptic integrals
  k2 =  4*a.*b./((a+b).^2 + c.^2);
  k = sqrt(k2);
 [K,E] = ellipke(k2); if(min(size(E))>=2);E=E(:,2);K=K(:,2); end;
  dEdk= (E./k-K./k);
  dKdk= (E./((1-k.^2).*k)-K./k);

  r=b;
  z=c;
  ak=k;

     g3=(a+r).^2+z.^2;
  dg3dc=2*z;

  dkdc=-(2*a.*r.*dg3dc)./(k.*g3.^2);
  
  F3=(2-k2).*K-2*E;
  dF3dc=(2-k2).*dKdk.*dkdc - 2*k.*dkdc.*K - 2*dEdk.*dkdc;

  dAdc= -2e-7*cur.*sqrt(a./r).*(dF3dc./k - F3.*dkdc./k2);

  clear ak k K E a b c r z g1 g2 g3 g1r g2r g3r akr Er Kr p1 p2 p2;

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

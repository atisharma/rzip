% A routine to create a filament plasma description based Tokamak model. Where the C-matrix is via flux and fields.
% Based on the RZIP concept and methods as outlined in Nuclear Fusion 39(5) 663. NO EIGEN MODE REDUCTION
%
% [A, B, C, D, curlyM, curlyR, plasma_outputs, IFLAG] = 
%           rzip(shot, time, Used_coils, plasma_inputs, Machine_name, force_make, extra_routine)
%
% Inputs
% ======
% shot,              Defines the plasma description. shot=0 implies no plasma. Requires "Machine_name_plasma.m"
% time,              If reading real machine data, the time at which to acquire the above.
% Used_coils,        Select which coils are relevant to the experiment eg [1 3 5], [1:18].
% plasma_inputs,     [Plasma resistance, d(Plasma resistance)/dr, d(beta_p)/dt, d(l_i)/dt]
% Machine_name,      A text-string refering to the local file "Machine_name" which describes the machine configuration. 
% force_make,        ignore workspace logic and always create files and save.
% extra_routine,     name of routine which can adjust curlyM and curlyR before output if required.
%
% Outputs 
% =======
% A, B, C, D         The model matrices. C is for Bpol probes and flux loops defined in the machine description.
% curlyM, curlyR     The full 'inductance' and 'resistance' matrices.
% plasma_outputs     [stabilising force, destabilising force, decay index, Shafranov vertical field,
%                     plasma averaged dbr/dz, vertical stability parameter alpha].
% IFLAG              A status flag for correct execution (0=fine)
%
% Notes
% -----
% Need Machine_name.m and Machine_name_plasma.m in the local directory. rzip home directory must be added by path.
% Full use of the code (plasma and plasma_less) produces two files. 
% ./Machine_name_plasma_less.mat which has the purely structural details within it.
% ./Machine_name_plasma_full.mat which has the coupling from the structure to the proposed plasma grid.
% Use lower case names for machines to ensure VMS compatibility. 
% Examine rzip.m for full details of what Machine_name.m and the plasma equilibrium needs to define.
% IMPORTANT  In Lausanne code had a tendency to run out of memory when creating a full model, usually during the plasma_full subroutine. 
%            Solution is to assume the plasma_less data is ok and rerun with force_make=99.
%
% Original version Dec 98
% updated to function call rzip and basic input/output list
% updated to expanded output/input and labelling convention changed
% updated to speed extensions, workspace checking etc
% updated to exact coupling from probes to plasma grid
% updated to use plasma specified on a non-retangular grid, implies use of A_plasma variable.
% updated to be asymmetric by default and to use logical names. Improved workspace logic at expense of data interpolation
% updated to use RIp as default 'radial' variable and curlyM symmetrised by using (radial equation)/IP
% updated to allow extra_routine and to force saving of files. Default of C matrix is for all states, not just coil currents.
%  
% John Wainwright Sep'99



% Defined by Machine_name.m
% =========================
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
% R_PS,		radial position
% Z_PS.		vertical position
% RES_PS.	resistance
% DR_PS,	width
% DZ_PS,	height
% LA_PS,	label to distinguish parts of structure: convention ??###F001
%                                                                     ^  ^  ^
%                                                                     |  |  L_ filament 1
%                                                                     |  L____ Filament 
%                                                                     L_______ name (5 chars 2/3 eg VV001)
% Diagnostics
% ----------- 
% R_BP,		radial position of Bpol probe
% Z_BP.		vertical position of Bpol probe
% TH_BP,        angle of Bpol probe
% R_FL,		radial position of flux loop
% Z_FL.		vertical position of flux loop
%
% Defined by Machine_name_plasma.m
% ================================
%
% plasma grid (grid over which plasma will be defined)
% ----------------------------------------------------
% R_plasma, 	  radial coordinate 
% Z_plasma,	  vertical coordinate 
% J_plasma,       plasma current density at point R_plasma, Z_plasma 
% A_plasma,       associated current area at point R_plasma, Z_plasma 
% BETA_P,         an external input for beta_p consistant with above
% LI,             an external input for li consistant with above


  function [A, B, C, D, curlyM, curlyR, plasma_outputs, IFLAG] = ...
               rzip_noeig(shot, time, Used_coils, plasma_inputs, Machine_name, force_make, extra_routine);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global variables associated with workspace 
  global Z_plasma_old R_plasma_old Machine_name_old Used_coils_old

% global variables associated with plasma grid
  global dbzdz_struct_grid  dbrdr_struct_grid  dbzdr_struct_grid    ...
           mut_struct_grid     br_struct_grid     bz_struct_grid    ...
            bprobe_to_grid       dbprobe_grid_dr    dbprobe_grid_dz ... 
             floop_to_grid        dfloop_grid_dr     dfloop_grid_dz              

% global variables associated with plasma less model
  global L_struct Res_struct Eigen_vectors ...
         R_struct Z_struct ...
         N_struct_fils N_active_coils N_active_fils N_passive_fils ...
         active_coil_ranges NT_AC floop_to_struct bprobe_to_struct

% flag for screen output on(1)/off(0)   
  screen_is_on = 1;

% Initialise situation dependent flags
% flag for decision on whether to make the plasma_less model (default no)
  make_plasma_less = 0;
% flag for decision on whether to make the machine-plasma grid links model (default no)
  make_plasma_full = 0;
% flag for decision on whether to use plasmma less workspace data, if it exists (default no)
  use_workspace_less = 0;
% flag for decision on whether to use plasmma full workspace data, if it exists (default no)
  use_workspace_full = 0;
% flag for decision on adding disturbance matrices (default no)
  make_disturb = 0;
% flag for whether a plasma is required (default no)
  plasma_less = 0; 
% flag for whether code has executed properly (0 is true)
  IFLAG = 0;

% Create two file name variables to the output files.
  plasma_less_file  = [Machine_name,'_plasma_less_noeig.mat'];
  plasma_full_file  = [Machine_name,'_plasma_full_noeig.mat'];

% if not a local file check if is in the path structure somewhere.
  remote_plasma_less_file=which(plasma_less_file);
  remote_plasma_full_file=which(plasma_full_file);
  if(length(remote_plasma_less_file)); plasma_less_file=remote_plasma_less_file;  end;
  if(length(remote_plasma_full_file)); plasma_full_file=remote_plasma_full_file;  end;

% default the initial 'old' variables
  if(~length(Machine_name_old));
   Machine_name_old = '   ';  Used_coils_old = [];  
  end

% if force_make is 99 then we need to perform the memory corrections
  if(force_make==99);   
   bz_struct_grid = [];   
   if(screen_is_on); disp('Memory error correction: assuming that plasma_less is done and ok'); end;
  end;

% Now we have implicit directory use, check file is present
  if(~length(feval('which',Machine_name))); 
   disp(['Cannot find machine description file ',Machine_name,'.m  . Are you in the correct directory?'])
   A=[];B=[];C=[];D=[];curlyM=[];curlyR=[];plasma_outputs=[];IFLAG=1; 
   return;
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE EXECUTION REQUIREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Is there a plasma (shot!=0)?
%=============================
 if(shot==0); 

% define the plasma_less flag
  plasma_less = 1; 
  if(screen_is_on); disp('Assuming plasmaless'); end; 

 else; 

% get the plasma equilibrium details (J_plasma, A_plasma, Z_plasma, R_plasma, BETA_P, LI, I_AC) 
% ---------------------------------------------------------------------------------------------
  [J_plasma, A_plasma, R_plasma, Z_plasma, BETA_P, LI, I_AC] = feval([Machine_name,'_plasma'],shot,time);

  plasma_resistance = plasma_inputs(1);   d_plasma_resistance_dr = plasma_inputs(2);   
  betap_dot         = plasma_inputs(3);   li_dot                 = plasma_inputs(4);

% if disturbances are defined must make disturbance matrices
  if(li_dot~=0 | betap_dot~=0) ; make_disturb=1 ;end;

% if we're using tcv then we should use all coils possible
  if(Machine_name(1:3)=='tcv');
   if(length(I_AC) < length(Used_coils)); Used_coils = 1:length(I_AC); end;
   if(max(Used_coils) > 19); Used_coils = 1:19; end;
   if(isnan(BETA_P)) ; A=[];B=[];C=[];D=[];curlyM=[];curlyR=[];plasma_outputs=[];IFLAG=1; return; end;
   if(isinf(BETA_P)) ; A=[];B=[];C=[];D=[];curlyM=[];curlyR=[];plasma_outputs=[];IFLAG=1; return; end;
  end;
 end;     % the plasma_less check loop

% If the machine has changed then need to ensure we load in the variables  
  if(any(Machine_name_old(1:3)~=Machine_name(1:3))); L_struct = []; bz_struct_grid = []; end;


% Are there plasma less variables in the workspace? (length(L_struct)>0) (NB don't use exist as initial global have size)
% will not reload unless i) Used_coils are different     
%======================================================================================================================== 
 if(length(L_struct) > 0);

% define the no loading flag
  use_workspace_less = 1;
 
% i) if the number of coils has changed, must reload (true for plasma and plasma_less cases)
% -------------------------------------------------------------------------------------------
   if(length(Used_coils_old)~=length(Used_coils));
    if(screen_is_on); disp('The number of coils has changed, must reload data'); end;
    use_workspace_less = 0; bz_struct_grid = []; 
   elseif(any(Used_coils_old ~= Used_coils));
    if(screen_is_on); disp('The number of coils has changed, must reload data'); end;
    use_workspace_less = 0; bz_struct_grid = []; 
   end;

 else;  % variables are not in the workspace but maybe we have them on file

% check whether the plasma-less model has been already created  
%------------------------------------------------------------- 
  fid = feval('fopen',plasma_less_file,'r');
  if(fid == -1); make_plasma_less = 1; else ; fclose(fid); end; 
  bz_struct_grid = []; % if we're reloading the structure we must reload the grid
 end;    % end the check for variables in the workspace  


% Are there plasma full variables in the workspace? (length(bz_struct_grid)>0) 
% will reload if plasma grids are different  
% -----------------------------------------------------------------------------
  if(length(bz_struct_grid) > 0);

% define the no loading flag
   use_workspace_full = 1;

  else % variables are not in the workspace but maybe we have them on file


% check whether the structure-grid data model has already been created 
%--------------------------------------------------------------------- 
   if(~plasma_less); 
    fid = feval('fopen',plasma_full_file,'r');
    if(fid==-1);
     make_plasma_full = 1; 
    else;
     fclose(fid);
     feval('load',plasma_full_file);
%     if(shot==1); load case1_used; 
%      if(length(R_plasma_old)~=length(R_plasma));
%       R_plasma_old=R_plasma_old(x); Z_plasma_old=Z_plasma_old(x);
%      end;
%     end;
    end;                  

%  i) if the plasma grids are different reload (only check plasma grid if there is a plasma)
% ------------------------------------------------------------------------------------------
    if(make_plasma_full~=1)  
     if(length(R_plasma_old)~=length(R_plasma) | length(Z_plasma_old)~=length(Z_plasma));
      make_plasma_full = 1; 
     elseif(any(R_plasma_old ~= R_plasma) | any(Z_plasma_old ~= Z_plasma));
      make_plasma_full = 1; 
     end;
    end;

   end;   % end of plasma_less check                  

  end;    % end the check for variables in the workspace  


% however all above is for naught if force_make is defined
  if(force_make);   
   if(force_make~=99); make_plasma_less = 1; use_workspace_less = 0; end;
   if(~plasma_less);   make_plasma_full = 1; use_workspace_full = 0; end; 
   if(screen_is_on);   disp('Forced file creation and saving');    end;
   force_text='Over-riding ';
  else
   force_text='Found no ';
  end;


% start clock
%------------
 tic 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there is no plasma less model we will have to create it %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if(make_plasma_less) 

  t_start_less = toc;
  f_start_less = flops;

% tell the user what's happening
%-------------------------------- 
  if(screen_is_on); disp([force_text,'saved plasmaless data; creating from scratch']); end;

% get the relevant machine details
%--------------------------------- 
  [R_AC,Z_AC,RES_AC,DR_AC,DZ_AC,LA_AC,NT_AC,R_PS,Z_PS,RES_PS,DR_PS,DZ_PS,LA_PS,R_BP,Z_BP,TH_BP,LA_BP,R_FL,Z_FL,LA_FL] = feval(Machine_name);

% To ensure consistancy work only with column inputs so make sure all arrays are columns;
%-----------------------------------------------------------------------------------------
   R_AC = check_is_column(R_AC);   Z_AC = check_is_column(Z_AC); RES_AC = check_is_column(RES_AC);
  DR_AC = check_is_column(DR_AC); DZ_AC = check_is_column(DZ_AC); LA_AC = check_is_column(LA_AC);  NT_AC = check_is_column(NT_AC);

   R_PS = check_is_column(R_PS);   Z_PS = check_is_column(Z_PS); RES_PS = check_is_column(RES_PS);
  DR_PS = check_is_column(DR_PS); DZ_PS = check_is_column(DZ_PS); LA_PS = check_is_column(LA_PS);

   R_BP = check_is_column(R_BP);   Z_BP = check_is_column(Z_BP);  TH_BP = check_is_column(TH_BP);  LA_BP = check_is_column(LA_BP);
   R_FL = check_is_column(R_FL);   Z_FL = check_is_column(Z_FL);  LA_FL = check_is_column(LA_FL);

% First of all we need to make decisions about how many actual coils there are as opposed to coil filaments. 
% This is information stored in LA (for label) initially we must work with filaments and combine to the coils.
%-------------------------------------------------------------------------------------------------------------

% active coils (N_active_coils = number of active coils, N_active_fils = number of active filaments)
%==============-------------------------------------------------------------------------------------
% seek out the different coils by name
 [N_active_coils,N_active_fils,active_coil_ranges,coil_names] = unique_name(LA_AC(:,1:3));

% passive coils (N_passive_coils = number of passive 'coils', N_passive_fils = number of passive filaments)
%===============-------------------------------------------------------------------------------------------
% seek out the different structures by name
 [N_passive_coils,N_passive_fils,passive_coil_ranges,passive_names] = unique_name(LA_PS(:,1:2));


 if(screen_is_on);
  disp(['Potentially ',num2str(N_active_coils),' Active coils and ',num2str(N_passive_coils),' Passive structure types']);
  fill=num2str(length(Used_coils));if(length(Used_coils)==N_active_coils);fill='all of the';end;
  disp(['will use ',fill,' Active coils and non-reduced structure']);
  disp(['There are ',num2str(length(R_BP)),' Bpol probes and ',num2str(length(R_FL)),' flux loops']);
 end;

% however will treat the separate passive structures as one structure (plotting use).
 N_passive_coils=1;passive_coil_ranges(1,1)=1;passive_coil_ranges(1,2)=N_passive_fils;

% so we can begin the inductance matrix calculation
%===================================================

% from the active filaments to the passive filaments
  for k=1:N_active_fils
   L_af_pf(k,1:N_passive_fils) = mutual_ind(R_AC(k),R_PS,Z_PS-Z_AC(k))';
  end;

% from the active filaments to the active filaments
  for k=1:N_active_fils
   L_af_af(k,1:N_active_fils) = mutual_ind(R_AC(k),R_AC,Z_AC-Z_AC(k))';
  end;

% from the passive filaments to the passive filaments
  for k=1:N_passive_fils
   L_pf_pf(k,1:N_passive_fils) = mutual_ind(R_PS(k),R_PS,Z_PS-Z_PS(k))';
  end;


% obviously there will be some difficulties we need to redefine the self inductances and to make sure that if coils are too 
% close we are accurate enough. (should be ok here). Also remember mutual_ind removes NaN and Inf and replaces with 0

%  L_af_af = L_af_af+diag(self_circ(R_AC,sqrt((DR_AC.*DZ_AC)/pi)));
  L_af_af = L_af_af+diag(self_sqre(R_AC,DR_AC,DZ_AC));
  L_pf_pf = L_pf_pf+diag(self_sqre(R_PS,DR_PS,DZ_PS));

% need to combine filaments to coils and include the number of turns  
%-------------------------------------------------------------------

% for the active/passive connections
  L_ac_ps = coil_from_filaments(L_af_pf.*[NT_AC*ones(1,N_passive_fils)],active_coil_ranges,N_active_coils,N_active_fils,0);

% for the active connections
% compact rows
  temp = coil_from_filaments(L_af_af.*[NT_AC*ones(1,N_active_fils)],active_coil_ranges,N_active_coils,N_active_fils,0)';
% then columns  
  L_ac_ac = coil_from_filaments(temp.*[NT_AC*ones(1,N_active_coils)],active_coil_ranges,N_active_coils,N_active_fils,0)';
 
% Memory saving
 clear L_af_af L_af_pf

 if(screen_is_on);disp(['Created the plasmaless L and R after ',num2str(toc-t_start_less),'s']);end;

% assume there is no compacting for the VV, if there is it would be achieved using a NT_PS variable

% Hijack the eigenmodes (see notes)
%=====================================

% define the do no reduction matrices
  Eigen_vectors = eye(length(L_pf_pf)); 
  Eigen_values  = L_pf_pf; 

% define the equivalent inductance matrix
%----------------------------------------
  L_struct = [        L_ac_ac                 L_ac_ps*Eigen_vectors ;
               [L_ac_ps*Eigen_vectors]'            Eigen_values    ];

% and the equivalent resistance matrix 
%--------------------------------------
  Res_ac_ac  = diag(coil_from_filaments(RES_AC, active_coil_ranges, N_active_coils, N_active_fils, 0));
  Res_struct = [          Res_ac_ac                          zeros(N_active_coils,N_passive_fils) ;
                 zeros(N_passive_fils,N_active_coils)         diag(RES_PS)   ]; 

% Memory saving
  clear L_pf_pf L_ac_ps inv_R inv_R_L Ee Es Ve V_o Eigen_values Vs Res_ac_ac 

% define vectors of filament positions (coils and structure)
  R_struct      = [R_AC;R_PS];
  Z_struct      = [Z_AC;Z_PS];
  N_struct_fils = length(R_struct);

% C Matrix information takes some construction (we're considering b-pol and flux loops)
%======================================================================================
 
  N_bprobe = length(R_BP); N_floop = length(R_FL); 

% Consider coupling of bpol probes (bp) to coils and structure (cs)
%------------------------------------------------------------------
  Bprcs=zeros(N_bprobe,N_struct_fils);Bpzcs=zeros(N_bprobe,N_struct_fils);
  for k=1:N_bprobe;
   Bprcs(k,:)  =   bfield_br(R_struct,R_BP(k),Z_BP(k)-Z_struct,1)'; 
   Bpzcs(k,:)  =   bfield_bz(R_struct,R_BP(k),Z_BP(k)-Z_struct,1)'; 
  end;

% convert to eigenmodes
  Bprcs = [Bprcs(:,1:N_active_fils) Bprcs(:,N_active_fils+1:N_active_fils+N_passive_fils)*Eigen_vectors];
  Bpzcs = [Bpzcs(:,1:N_active_fils) Bpzcs(:,N_active_fils+1:N_active_fils+N_passive_fils)*Eigen_vectors];
 
% but have to create the net effects of the links to the coils

  Bprcs  = Bprcs'.*[[NT_AC;ones(N_passive_fils,1)]*ones(1,N_bprobe)];
  Bpzcs  = Bpzcs'.*[[NT_AC;ones(N_passive_fils,1)]*ones(1,N_bprobe)]; 
  bprcs   = coil_from_filaments(Bprcs,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils)';
  bpzcs   = coil_from_filaments(Bpzcs,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils)';

% Create total field(bpt) and angle(bpa)
  bptcs   = sqrt(bprcs.*bprcs+bpzcs.*bpzcs); 
  bpacs   = atan2(bpzcs,bprcs); 

% compare with angle of coil
  net_angle = bpacs - [TH_BP*ones(1,N_active_coils+N_passive_fils)];

% so the final effect from coils/structure to bpols is
  bprobe_to_struct = cos(net_angle).*bptcs;


% and now the flux probes (fl) (which is simply via mutual inductance F=LI)
%---------------------------------------------------------------------------
% first to the coils and structure (cs)
  Flcs = zeros(N_floop,N_struct_fils);
  for k=1:N_floop;
   Flcs(k,:) = mutual_ind(R_struct,R_FL(k),Z_FL(k)-Z_struct)'; 
  end;
% convert to eigenmodes
  Flcs = [Flcs(:,1:N_active_fils) Flcs(:,N_active_fils+1:N_active_fils+N_passive_fils)*Eigen_vectors];
  
% but have to create the net effects of the links to the coils
  Flcs              = Flcs'.*[[NT_AC;ones(N_passive_fils,1)]*ones(1,N_floop)];
  floop_to_struct   = coil_from_filaments(Flcs,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils)';

% Save the necessary outputs for later
  old_clock=clock;
  feval('save',plasma_less_file, ...
               'N_passive_fils','N_passive_fils','N_active_coils','N_active_fils', ...  
               'active_coil_ranges','NT_AC','R_struct','Z_struct','N_struct_fils', ...
               'bprobe_to_struct','floop_to_struct','L_struct','Res_struct', ...
               'Eigen_vectors','old_clock');
  if(screen_is_on);
   disp(['Creating the plasmaless model required ',num2str(toc-t_start_less),'s and ',num2str(flops-f_start_less,8),' flops']);
  end;

  clear Flcs bpacs bptcs bprcs bpzcs net_angle Bprcs Bpzcs norm_mat

 else;   % read in the old plasma_less data
         %---------------------------------

% tell the user what's happening
%-------------------------------- 
  if(use_workspace_less)
   if(screen_is_on);disp('Using workspace plasma less data');end;
  else
   feval('load',plasma_less_file);
   if(screen_is_on);
    creation_date=[num2str(old_clock(3)),'/',num2str(old_clock(2)),'/',num2str(old_clock(1))];
    if(old_clock(4)<10);fil1='0';else;fil1='';end;
    if(old_clock(5)<10);fil2=':0';else;fil2=':';end;
    creation_time=[fil1,num2str(old_clock(4)),fil2,num2str(old_clock(5))];
    disp(['Using saved plasma_less data from ',creation_time,' on the ',creation_date])
   end;
  end;

 end;   % make_plasma_less loop
        %----------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the appropriate coils %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the number of used coils (N_Used_coils) and define the matrix of used circuits
% --------------------------------------------------------------------------------------
 if(use_workspace_less)
                 N_Used_coils = length(Used_coils);
                Used_circuits = [1:(N_Used_coils+N_passive_fils)];
 else
                Used_circuits = [Used_coils,N_active_coils+1:N_active_coils+N_passive_fils];
                 N_Used_coils = length(Used_coils);
                     L_struct = L_struct(Used_circuits,Used_circuits);
                   Res_struct = Res_struct(Used_circuits,Used_circuits);
             bprobe_to_struct = bprobe_to_struct(:,Used_circuits);
              floop_to_struct = floop_to_struct(:,Used_circuits);
 end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there are no links from the machine to the grid the will have to create them %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if(make_plasma_full) 

  t_start_full = toc;
  f_start_full = flops;
 N_grid_points = length(R_plasma); 

% tell the user what's happening
%-------------------------------- 
 if(screen_is_on);    disp([force_text,'saved plasmafull data; creating from scratch']);     end;

% however at this stage if we haven't performed the plasma_less
% model then we need to reobtain the machine data
%------------------------------------------------
 if(~make_plasma_less) ;
  [R_AC,Z_AC,RES_AC,DR_AC,DZ_AC,LA_AC,NT_AC,R_PS,Z_PS,RES_PS,DR_PS,DZ_PS,LA_PS,R_BP,Z_BP,TH_BP,LA_BP,R_FL,Z_FL,LA_FL] = feval(Machine_name);
  N_bprobe = length(R_BP); N_floop = length(R_FL); 
  R_struct=[R_AC;R_PS];  Z_struct=[Z_AC;Z_PS];  N_struct_fils=length(R_struct);
 end;

% we need to define various interactions of the structure to the plasma grid. Assume the number of grid points is the larger
% and therefore optimise the calculations to this.  We need Bz, Br, Mut(a/p,g) dBzdr dBrdr dBzdz (the bulk of the code works in this bit) 

   TMP = zeros(N_struct_fils,N_grid_points);
   for k=1:N_struct_fils;
    TMP(k,:) = bfield_bz(R_struct(k),R_plasma,Z_plasma-Z_struct(k),1)'; 
   end;

% with an eigenmode analysis we need to adjust these fields to be scaled by the eigenvectors. 
% Ie we convert [X_af X_ps]*[I_af I_ps] to  [X]*[I_af I_e]  X=[X_af X_ps*E] 
% (bear in mind the transposition  due to definition of eg bz to grid) 
   TMP = [TMP(1:N_active_fils,:);[TMP(N_active_fils+1:N_active_fils+N_passive_fils,:)'*Eigen_vectors]'];

% but have to create the net effects of the links to the coils. Must scale by number of turns
   TMP = [[NT_AC;ones(N_passive_fils,1)]*ones(1,N_grid_points)].*TMP;
   bz_struct_grid = coil_from_filaments(TMP,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils);
   if(screen_is_on);disp('Finished with bz_struct_grid');end;

   TMP = zeros(N_struct_fils,N_grid_points);
   for k=1:N_struct_fils;
    TMP(k,:) = bfield_br(R_struct(k),R_plasma,Z_plasma-Z_struct(k),1)'; 
   end;
   TMP = [TMP(1:N_active_fils,:);[TMP(N_active_fils+1:N_active_fils+N_passive_fils,:)'*Eigen_vectors]'];
   TMP = [[NT_AC;ones(N_passive_fils,1)]*ones(1,N_grid_points)].*TMP;
   br_struct_grid    = coil_from_filaments(TMP,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils);
   if(screen_is_on);disp('Finished with br_struct_grid');end;

   TMP = zeros(N_struct_fils,N_grid_points);
   for k=1:N_struct_fils;
    TMP(k,:) = mutual_ind(R_struct(k),R_plasma,Z_plasma-Z_struct(k))'; 
   end;
   TMP = [TMP(1:N_active_fils,:);[TMP(N_active_fils+1:N_active_fils+N_passive_fils,:)'*Eigen_vectors]'];
   TMP = [[NT_AC;ones(N_passive_fils,1)]*ones(1,N_grid_points)].*TMP;
   mut_struct_grid    = coil_from_filaments(TMP,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils);
   if(screen_is_on);disp('Finished with mut_struct_grid');end;

   TMP = zeros(N_struct_fils,N_grid_points);
   for k=1:N_struct_fils;
    TMP(k,:) = bfield_dbrdr(R_struct(k),R_plasma,Z_plasma-Z_struct(k),1)'; 
   end;
   TMP = [TMP(1:N_active_fils,:);[TMP(N_active_fils+1:N_active_fils+N_passive_fils,:)'*Eigen_vectors]'];
   TMP = [[NT_AC;ones(N_passive_fils,1)]*ones(1,N_grid_points)].*TMP;
   dbrdr_struct_grid    = coil_from_filaments(TMP,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils);
   if(screen_is_on);disp('Finished with dbrdr_struct_grid');end;

   TMP = zeros(N_struct_fils,N_grid_points);
   for k=1:N_struct_fils;
    TMP(k,:) = bfield_dbzdz(R_struct(k),R_plasma,Z_plasma-Z_struct(k),1)'; 
   end;
   TMP = [TMP(1:N_active_fils,:);[TMP(N_active_fils+1:N_active_fils+N_passive_fils,:)'*Eigen_vectors]'];
   TMP = [[NT_AC;ones(N_passive_fils,1)]*ones(1,N_grid_points)].*TMP;
   dbzdz_struct_grid    = coil_from_filaments(TMP,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils);
   if(screen_is_on);disp('Finished with dbzdz_struct_grid');end;

   TMP = zeros(N_struct_fils,N_grid_points);
   for k=1:N_struct_fils;
    TMP(k,:) = bfield_dbzdr(R_struct(k),R_plasma,Z_plasma-Z_struct(k),1)'; 
   end;
   TMP = [TMP(1:N_active_fils,:);[TMP(N_active_fils+1:N_active_fils+N_passive_fils,:)'*Eigen_vectors]'];
   TMP = [[NT_AC;ones(N_passive_fils,1)]*ones(1,N_grid_points)].*TMP;
   dbzdr_struct_grid = coil_from_filaments(TMP,active_coil_ranges,N_active_coils,N_active_fils,N_passive_fils);
   clear TMP Eigen_vectors
   if(screen_is_on);
    disp('Finished with dbzdr_struct_grid');
    disp(['Creating the main matrices took ',num2str(toc-t_start_full),' seconds']);
   end;

% plasma grid mapping to the diagnostics. Using exact fields and using new routines
%===================================================================================

% first get coupling of bpol probes to plasma grid (g)
%-----------------------------------------------------
  dbprobe_grid_dr = zeros(N_bprobe,N_grid_points); dbprobe_grid_dz = zeros(N_bprobe,N_grid_points);
          bpdrdrg = zeros(N_bprobe,N_grid_points);         bpdzdrg =  zeros(N_bprobe,N_grid_points);
          bpdzdzg = zeros(N_bprobe,N_grid_points);

  for k=1:N_bprobe;
   dbprobe_grid_dr(k,:)  =    bfield_br(R_plasma,R_BP(k),Z_BP(k)-Z_plasma,1)'; 
   dbprobe_grid_dz(k,:)  =    bfield_bz(R_plasma,R_BP(k),Z_BP(k)-Z_plasma,1)'; 
           bpdzdrg(k,:)  =    bfield_dbzda(R_plasma,R_BP(k),Z_BP(k)-Z_plasma,1)'; 
           bpdzdzg(k,:)  =    bfield_dbzdc(R_plasma,R_BP(k),Z_BP(k)-Z_plasma,1)'; 
           bpdrdrg(k,:)  =    bfield_dbrda(R_plasma,R_BP(k),Z_BP(k)-Z_plasma,1)'; 
           bpdrdzg(k,:)  =    bfield_dbrdc(R_plasma,R_BP(k),Z_BP(k)-Z_plasma,1)'; 
  end;

% Create total field(bpt) and angle(bpa)
  bptg   = sqrt(dbprobe_grid_dr.*dbprobe_grid_dr+dbprobe_grid_dz.*dbprobe_grid_dz); 
  bpag   = atan2(dbprobe_grid_dz,dbprobe_grid_dr); 

% create the linear terms (see the algebra)
  dthetadr = (dbprobe_grid_dr.*bpdzdrg-dbprobe_grid_dz.*bpdrdrg)./bptg;
  dthetadz = (dbprobe_grid_dr.*bpdzdzg-dbprobe_grid_dz.*bpdrdzg)./bptg;
  dmagdr   = (dbprobe_grid_dr.*bpdrdrg+dbprobe_grid_dz.*bpdzdrg)./bptg;
  dmagdz   = (dbprobe_grid_dr.*bpdrdzg+dbprobe_grid_dz.*bpdzdzg)./bptg;

% compare with angle of coil
  net_angle = [TH_BP*ones(1,N_grid_points)] - bpag;

% so the final effect from grid to bpols is
   bprobe_to_grid = cos(net_angle).*bptg;
  dbprobe_grid_dr = sin(net_angle).*dthetadr+cos(net_angle).*dmagdr;
  dbprobe_grid_dz = sin(net_angle).*dthetadz+cos(net_angle).*dmagdz;
 
% and now the flux probes (fl) (which is simply via mutual inductance F=LI)
%--------------------------------------------------------------------------
  floop_to_grid = zeros(N_floop,N_grid_points);  VPr=zeros(N_floop,N_grid_points); VPz=zeros(N_floop,N_grid_points); 

  for k=1:N_floop;
   floop_to_grid(k,:)  =     mutual_ind(R_plasma,R_FL(k),Z_FL(k)-Z_plasma)'; 
             VPr(k,:)  =     vecpot_dada(R_plasma,R_FL(k),Z_FL(k)-Z_plasma,1)'; 
             VPz(k,:)  =     vecpot_dadc(R_plasma,R_FL(k),Z_FL(k)-Z_plasma,1)'; 
  end;

  dfloop_grid_dr  =  2*pi*[R_FL*ones(1,N_grid_points)].*VPr;
  dfloop_grid_dz  =  2*pi*[R_FL*ones(1,N_grid_points)].*VPz;

  R_plasma_old = R_plasma; Z_plasma_old = Z_plasma; old_clock = clock;

  feval('save', plasma_full_file, ...
                 'floop_to_grid','dfloop_grid_dr','dfloop_grid_dz', ...
                 'bprobe_to_grid','dbprobe_grid_dr','dbprobe_grid_dz', ... 
                 'bz_struct_grid','br_struct_grid', ...
                 'dbzdr_struct_grid','dbzdz_struct_grid','dbrdr_struct_grid', ... 
                 'mut_struct_grid','old_clock','R_plasma_old','Z_plasma_old');

  if(screen_is_on);
   disp(['Creating the links to the plasma grid required ',num2str(toc-t_start_full),'s and ',num2str(flops-f_start_full,8),' flops']);
  end;

 else;   % read in the old  plasma_full data

% tell the user what's happening
%-------------------------------- 
  if(use_workspace_full)
   if(screen_is_on);disp('Using workspace plasma full data');end;
  else
   if(~plasma_less)
    feval('load',plasma_full_file);

%    if(shot==1); load case1_used; 
%     R_plasma_old=R_plasma_old(x);
%     dbprobe_grid_dz=dbprobe_grid_dz(:,x);
%     floop_to_grid= floop_to_grid(:,x);      
%     Z_plasma_old=Z_plasma_old(x);
%     dbrdr_struct_grid=dbrdr_struct_grid(:,x);
%     mut_struct_grid=mut_struct_grid(:,x);     
%     bprobe_to_grid=bprobe_to_grid(:,x);
%     dbzdr_struct_grid=dbzdr_struct_grid(:,x);
%     br_struct_grid=br_struct_grid(:,x);
%     dbzdz_struct_grid=dbzdz_struct_grid(:,x);   
%     bz_struct_grid=bz_struct_grid(:,x);
%     dfloop_grid_dr=dfloop_grid_dr(:,x);       
%     dbprobe_grid_dr=dbprobe_grid_dr(:,x);
%     dfloop_grid_dz=dfloop_grid_dz(:,x);      
%    end;

    if(screen_is_on);
     creation_date=[num2str(old_clock(3)),'/',num2str(old_clock(2)),'/',num2str(old_clock(1))];
     if(old_clock(4)<10);fil1='0';else;fil1='';end;
     if(old_clock(5)<10);fil2=':0';else;fil2=':';end;
     creation_time=[fil1,num2str(old_clock(4)),fil2,num2str(old_clock(5))];
     disp(['Using saved plasma_full data from ',creation_time,' on the ',creation_date])
    end;
   end;
  end;

 end;   % make_plasma_full


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the appropriate coils %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if(~plasma_less);
      bz_struct_grid =     bz_struct_grid(Used_circuits,:);
      br_struct_grid =     br_struct_grid(Used_circuits,:);
   dbzdr_struct_grid =  dbzdr_struct_grid(Used_circuits,:);
   dbzdz_struct_grid =  dbzdz_struct_grid(Used_circuits,:);
   dbrdr_struct_grid =  dbrdr_struct_grid(Used_circuits,:);
     mut_struct_grid =    mut_struct_grid(Used_circuits,:);
 end;


% set up the 'old' variables for use in comparing present and previous execution
  Machine_name_old = Machine_name;  Used_coils_old =  Used_coils;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make rzip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% And then make a proper plasmaless or rzip model
%================================================

 if(plasma_less)

  curlyM = L_struct ;
  curlyR = Res_struct ;
   Cbpol = bprobe_to_struct;
    Cpsi = floop_to_struct;
   Cipol = [eye(N_Used_coils) zeros(N_Used_coils,N_passive_fils)];
   Cstat = [eye(N_Used_coils+N_passive_fils)];

% set the plasma output terms to nonexistent

  plasma_outputs=[];

 else; %%%%%%%%%%%%% else of plasma less check

% define mu0
  mu0 = 4e-7*pi;

% define the equilibrium circuit currents
 struct_currents = [I_AC(Used_coils);zeros(N_passive_fils,1)];

% plasma current ~ int J_plasma A_plasma
 IP    = sum(J_plasma.*A_plasma);

% need to have read in Beta_p  externally
 betap = BETA_P;
    li = LI;

% plasma mutual inductance gradients
%-----------------------------------
 dL_plas_struct_dr  =  2*pi*(bz_struct_grid*[R_plasma.*J_plasma.*A_plasma])/IP;
 dL_plas_struct_dz  = -2*pi*(br_struct_grid*[R_plasma.*J_plasma.*A_plasma])/IP;

% plasma field gradients
%----------------------- 
 dbzdr_eqm    = struct_currents'*dbzdr_struct_grid;
 dbzdr_plas   = (dbzdr_eqm*[J_plasma.*A_plasma])/IP;
 dbrdz_plas   =  dbzdr_plas;
 rdbrdz_plas  = (dbzdr_eqm*[R_plasma.*J_plasma.*A_plasma])/IP;

% Useful equilibrium terms
 alpha      = 2*pi*rdbrdz_plas/IP;
 bz_eqm     = bz_struct_grid'*struct_currents; 
 br_eqm     = br_struct_grid'*struct_currents; 
 bz_plas    = bz_eqm'*[J_plasma.*A_plasma]/IP;
 br_plas    = br_eqm'*[J_plasma.*A_plasma]/IP;
 rbz_plas   = bz_eqm'*[R_plasma.*J_plasma.*A_plasma]/IP;
 rbr_plas   = br_eqm'*[R_plasma.*J_plasma.*A_plasma]/IP;
 Decay_plas = -dbzdr_eqm*(R_plasma.*J_plasma.*A_plasma./bz_eqm)/IP;
 Fd         = -2*pi*IP*rdbrdz_plas;
 index_ps   = N_Used_coils+1:N_Used_coils+N_passive_fils;
 Fs         = dL_plas_struct_dz(index_ps)'*inv(L_struct(index_ps,index_ps))*dL_plas_struct_dz(index_ps)*IP*IP; 

 plasma_outputs = [Fs, Fd, Decay_plas, bz_plas, dbrdz_plas, alpha];

% for destabilising force analysis
 save destab_noeig.mat dL_plas_struct_dz L_struct IP index_ps


% plasma properties
%------------------  
 R0       =  sum(R_plasma'*[J_plasma.*A_plasma])/IP; 
 gamma0   = -4*pi*R0*bz_plas/mu0/IP; 
 fo       = gamma0-betap-1/2;
 Lp0      = mu0*R0*fo;
 mut_plas = mut_struct_grid*[J_plasma.*A_plasma]/IP;

% asymmetric coefficiients
 dbrdr_eqm    = struct_currents'*dbrdr_struct_grid;
 rdbrdr_plas  = (dbrdr_eqm*[R_plasma.*J_plasma.*A_plasma])/IP;
 dbzdz_eqm    = struct_currents'*dbzdz_struct_grid;
 rdbzdz_plas  = (dbzdz_eqm*[R_plasma.*J_plasma.*A_plasma])/IP;

% display useful coefficients
%----------------------------
 if(screen_is_on);disp(['bz_plas,alpha = ', num2str(bz_plas),'  ',num2str(alpha)]);end;


% create the model 'inductance' and 'resistance' matrices (using R Ip as 'radial' term and 'radial' equation divided by IP to symmetrise curlyM)
%===============================================================================================================================================
 M11 = L_struct          ;M12 = dL_plas_struct_dz   ;M13 = dL_plas_struct_dr                            ;M14 = mut_plas;
 M21 = dL_plas_struct_dz';M22 = -2*pi*rdbrdz_plas/IP;M23 = -2*pi*rdbrdr_plas/IP                         ;M24 = 0;
 M31 = dL_plas_struct_dr';M32 =  2*pi*rdbzdz_plas/IP;M33 = (mu0*IP/(2*R0)+2*pi*(rdbrdz_plas+bz_plas))/IP;M34 = (mu0*IP*gamma0+2*pi*rbz_plas)/IP;
 M41 = mut_plas'         ;M42 = 0                   ;M43 = mu0*(1+fo)+2*pi*rbz_plas/IP                  ;M44 = Lp0;
 

 curlyM = [ M11  M12  M13 M14; 
            M21  M22  M23 M24; 
            M31  M32  M33 M34
            M41  M42  M43 M44];

 curlyR = [ Res_struct  zeros(length(Res_struct),3);
            zeros(2,length(curlyM));
            zeros(1,length(curlyM)-2) d_plasma_resistance_dr/IP plasma_resistance];   


% rigid plasma movement contribution to bpol probes 
  bp_plas =   bprobe_to_grid*[J_plasma.*A_plasma]/IP;
  dbpdr   =  dbprobe_grid_dr*[J_plasma.*A_plasma]/IP;
  dbpdz   =  dbprobe_grid_dz*[J_plasma.*A_plasma]/IP;
% rigid plasma movement contribution to flux loops 
  fl_plas =   floop_to_grid*[J_plasma.*A_plasma]/IP;
  dfldr   =  dfloop_grid_dr*[J_plasma.*A_plasma]/IP;
  dfldz   =  dfloop_grid_dz*[J_plasma.*A_plasma]/IP;

  Cr    = [zeros(1,N_Used_coils+N_passive_fils) 0 1 0];
  Cip   = [zeros(1,N_Used_coils+N_passive_fils) 0 0 1];
  Czip  = [zeros(1,N_Used_coils+N_passive_fils) 1 0 0];
  Cipol = [eye(N_Used_coils) zeros(N_Used_coils,N_passive_fils+3)];
  Cstat = [eye(N_Used_coils+N_passive_fils+3)];
  Cbpol = [bprobe_to_struct dbpdz dbpdr bp_plas];
  Cpsi  = [floop_to_struct  dfldz dfldr fl_plas];

 end; %%%%%%%%%%%%% end of plasma less check


% call extra manipulation routine if defined. Allows opportunity to manipulate curlyM and curlyR before use.
 if(length(extra_routine) > 0); 
%  feval(extra_routine) ;
  eval([extra_routine,';']);
 end


% Make the A B C and D matrices
% -----------------------------
 minv = inv(curlyM);
 extra_zero=0; if(~plasma_less); extra_zero=3; end;
 volt = [diag(ones(N_Used_coils,1)); zeros(N_passive_fils+extra_zero,N_Used_coils)]; 

 A =-minv*curlyR;
 B = minv*volt;
 C = [Cbpol; Cpsi; Cstat];
 D = zeros(size(C,1),size(B,2));

% If we have disturbances then we need to rearrange things
 if(make_disturb)
  e_dist = -mu0*IP*IP/2;
  E      = minv*[zeros(N_Used_coils+N_passive_fils,2); e_dist e_dist/2; 0 0];
% Rearange states so that Z=x-Ew ; and V=[v_coils beta_dot li_dot];
  B      = [B A*E];
  D      = [D C*E];
 end;

 if(screen_is_on);
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

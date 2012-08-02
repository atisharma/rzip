 
 if(~length(findstr('rzip',path))) ; addpath('/home/ps/sunfs1/d4/jw1/Matlab/rzip'); end;
 if(~length(findstr('rzip/iter',path))) ; addpath('/home/ps/sunfs1/d4/jw1/Matlab/rzip/iter'); end;
 if(~length(findstr('rzip/subroutines',path))) ; addpath('/home/ps/sunfs1/d4/jw1/Matlab/rzip/subroutines'); end;

 [R_AC,Z_AC,RES_AC,DR_AC,DZ_AC,LA_AC,NT_AC,R_PS,Z_PS,RES_PS,DR_PS,DZ_PS,LA_PS,R_BP,Z_BP,TH_BP,LA_BP,R_FL,Z_FL,LA_FL] = tcv;

% define a gross list of alternative names (LA_AC(:,1:3))
 [N_active_coils,N_active_fils,active_coil_ranges,coil_names] = unique_name(LA_AC(:,1:3));

% identify coils in a group
  coil_nums = str2num(LA_AC(:,4:5));

% however a 'coil' could be a compound of coils as well as filaments.
 for k = 1:N_active_coils
  coil_range = active_coil_ranges(k,1):active_coil_ranges(k,2);
  coils_in_coil(k,:) = max(coil_nums(coil_range));
 end;

% loop over coils and defilament them
  nco = 1;
 for k = 1:N_active_coils-1
  coil_range=active_coil_ranges(k,1):active_coil_ranges(k,2);
  for l = 1:coils_in_coil(k)
   coil_sub = find(coil_nums(coil_range)==l);
   N_fil_co = length(coil_sub);
   R_MAX = max(R_AC(coil_range(coil_sub))+DR_AC(coil_range(coil_sub))/2);
   R_MIN = min(R_AC(coil_range(coil_sub))-DR_AC(coil_range(coil_sub))/2);
   Z_MAX = max(Z_AC(coil_range(coil_sub))+DZ_AC(coil_range(coil_sub))/2);
   Z_MIN = min(Z_AC(coil_range(coil_sub))-DZ_AC(coil_range(coil_sub))/2);
   BIG_AREA = abs((R_MAX-R_MIN)*(Z_MAX-Z_MIN));
   BIG_RADI = (R_MAX+R_MIN)/2;
   BIG_RESI = 2*pi*BIG_RADI/BIG_AREA;
   FIL_AREA = DZ_AC(coil_range(coil_sub))'*DR_AC(coil_range(coil_sub));
   FIL_RES  = 2*pi*R_AC(coil_range(coil_sub))./( DR_AC(coil_range(coil_sub)).*DZ_AC(coil_range(coil_sub)) ) ; 
   FIL_RESI = sum( 2*pi*R_AC(coil_range(coil_sub))./( DR_AC(coil_range(coil_sub)).*DZ_AC(coil_range(coil_sub)) ) ); 
   BA(nco,1:2)=[BIG_AREA FIL_AREA];
   BR(nco,1:3)=[BIG_RESI FIL_RESI 1/sum(1./FIL_RES)];
   disp(['Areas and differences ',num2str(BIG_AREA),' ',num2str(FIL_AREA),' ',num2str(BIG_AREA/FIL_AREA)])
   nco = nco+1;
  end;
 end;
 






l=001:144;sum(RES_AC(l).*NT_AC(l).^2)/sum(2*pi*RES_TCV(1:16,2)./RES_TCV(1:16,3))/9   
l=145:176;sum(RES_AC(l).*NT_AC(l).^2)/sum(2*pi*RES_TCV(17,2)./RES_TCV(17,3))/32
l=209:220;sum(RES_AC(l).*NT_AC(l).^2)/sum(2*pi*RES_TCV(19,2)./RES_TCV(19,3))/12
l=233:240;sum(RES_AC(l).*NT_AC(l).^2)/sum(2*pi*RES_TCV(21,2)./RES_TCV(21,3))/8  
l=249:284;sum(RES_AC(l).*NT_AC(l).^2)/sum(2*pi*RES_TCV(23,2)./RES_TCV(23,3))/36
l=537:572;sum(RES_AC(l).*NT_AC(l).^2)/sum(2*pi*RES_TCV(31,2)./RES_TCV(31,3))/36 
l=825:825;sum(RES_AC(l).*NT_AC(l).^2)/sum(2*pi*RES_TCV(39,2)./RES_TCV(39,3))/1    
  

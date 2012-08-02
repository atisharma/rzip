% a script to print out the lines relevant to getting PROTEUS
% working for an equilibrium

  function for_pro(shot)

% first load in the appropriate data
cd data
eval(['load shot_',num2str(shot)]);

% define the plasma current
  IP = (J_PL'*A_PL);

% current set up is for p' = 1-x and f' = (1-x) + cf_rel*(1-x)^2
  if(length(find(TTC))~=2 & length(find(PPC))~=1)
   disp('Coefficients odd')
   [TTC' ; PPC']
   break;
  else;
   cp = PPC(2)/IP;
   cf = TTC(2)/IP;
   cf_rel = TTC(3)/TTC(2);
  end;

  disp(' CHANGES in fort.14')
  disp(' ')
  disp(['    plc = ',num2str(IP),', cp = ',num2str(cp),', cf = ',num2str(cf),','])
  disp(['    rzaxis(1)=',num2str((R_PL.*J_PL)'*A_PL/IP),',  rzaxis(2)= ',num2str((Z_PL.*J_PL)'*A_PL/IP),','])

  fJ = find(J_PL);
  [min_R,I_R] = min(R_PL(fJ));

  disp(['    rzboun(1)=',num2str(R_PL(fJ(I_R))),',  rzboun(2)= ',num2str(Z_PL(fJ(I_R))),','])
  disp(' ')

% try changing currents to compensate for the number of turns things
   nt=[144/143;1;ones(8,1)*36/34.1;ones(8,1)*36/35.9;.99;ones(5,1)];
  for k=1:length(I_AC)
   disp(['  amp(',num2str(k),') = ',num2str(I_AC(k)*nt(k)),','])
%   disp(['  amp(',num2str(k),') = ',num2str(I_AC(k)),','])
  end

  disp(' ')
  disp(' ')
  disp('Symmetric amps if preferred')
  disp(' ')


% symmetric amps
  CPS = [3 10;4 9;5 8;6 7;11 18;12 17;13 16; 14 15];
  amp =I_AC;
  for k=1:8
   a_sym = mean(I_AC(CPS(k,:)));
   amp(CPS(k,1))=a_sym;
   amp(CPS(k,2))=a_sym;
  end
  for k=1:length(I_AC)
   disp(['  amp(',num2str(k),') = ',num2str(amp(k)),','])
  end



  disp(' ')
  disp(' ')
  disp(' CHANGES in user.f')
  disp(' ')
  
  disp(['            cf_rel     = ',num2str(cf_rel)])

  disp(' ')
  disp(' ')
  disp(' expected fluxes axis|boundary')


  disp(['   ',num2str(max(max(PSI_PL+PSI_0))/2/pi)])
  disp(['   ',num2str(PSI_0/2/pi)])
  
 cd ..  
  
  return;
 

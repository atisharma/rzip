% Check we're not too far off linearised proteus with rzip
 
 if(~length(findstr('rzip',path))) ; addpath('/home/ps/sunfs1/d4/jw1/Matlab/rzip'); end;
 if(~length(findstr('rzip/iter',path))) ; addpath('/home/ps/sunfs1/d4/jw1/Matlab/rzip/iter'); end;

  used_coils = [1:19]; 
  nuc        = length(used_coils);
  nue        = 38;
  nus        = nuc+nue;

% get an rzip model
 [A, B, C, D, curlyM, curlyR, plasma_outputs, IFLAG] = rzip_v3(100,1,used_coils,nue,[0 0 0 0],'tcv',0,'');

% eliminate plasma variables to get a comparison
  Lrzip = curlyM(1:nus,1:nus) - curlyM(1:nus,nus+1:nus+3)*inv(curlyM(nus+1:nus+3,nus+1:nus+3))*curlyM(nus+1:nus+3,1:nus);
  Rrzip = curlyR(1:nus,1:nus); 

  Arzip = -inv(Lrzip)*Rrzip;
  Brzip = inv(Lrzip);
  Brzip = Brzip(:,used_coils);
  Crzip = -inv(curlyM(nus+1:nus+3,nus+1:nus+3))*curlyM(nus+1:nus+3,1:nus); Crzip = [Crzip(2,:);Crzip(1,:)];

% eliminate with fixed Ip
  Lrzif = curlyM(1:nus,1:nus) - curlyM(1:nus,nus+1:nus+2)*inv(curlyM(nus+1:nus+2,nus+1:nus+2))*curlyM(nus+1:nus+2,1:nus);
  Rrzif = curlyR(1:nus,1:nus); 

  Arzif = -inv(Lrzif)*Rrzif;
  Brzif = inv(Lrzif);
  Brzif = Brzif(:,used_coils);
  Crzif = -inv(curlyM(nus+1:nus+2,nus+1:nus+2))*curlyM(nus+1:nus+2,1:nus);
  Crzif = [Crzif(2,:);Crzif(1,:)];


% get a linearised PROTEUS
 ncircs = 77;

 Lp = zeros(ncircs,ncircs);
 Lm = zeros(ncircs,ncircs);
 Cp = zeros(2,ncircs);
 Cm = zeros(2,ncircs);

 filepathp  = ['/scratch/jw1/TCV/Fluxp/'];
 filepathm  = ['/scratch/jw1/TCV/Fluxm/'];
 for k = 1:ncircs
  filenamec = ['C',num2str(k)];
  filenamel = ['L',num2str(k)];

  fid = fopen([filepathp,filenamel],'r');
  for l = 1:ncircs
   linein = fgetl(fid);  numline = str2num(linein);
   Lp(k,l) = numline(3);
  end;
  fclose(fid);

  fid = fopen([filepathm,filenamel],'r');
  for l = 1:ncircs
   linein = fgetl(fid);  numline = str2num(linein);
   Lm(k,l) = numline(3);
  end;
  fclose(fid);

  fid = fopen([filepathp,filenamec],'r');
  linein = fgetl(fid);  numline = str2num(linein);
  Cp(1,k) = numline(3);
  Cp(2,k) = numline(4);
  fclose(fid);

  fid = fopen([filepathm,filenamec],'r');
  linein = fgetl(fid);  numline = str2num(linein);
  Cm(1,k) = numline(3);
  Cm(2,k) = numline(4);
  fclose(fid);

 end;

 Lstar = 0.5*(Lp(1:ncircs,1:ncircs)+Lm(1:ncircs,1:ncircs));
 
 Rstar = zeros(ncircs,ncircs);
 load /scratch/jw1/TCV/RES_TCV
 used_els = [1 16; 17 22; [23:38]' [23:38]';39 44;[45:102]' [45:102]'];
 for k=1:length(used_els);
  Rstar(k,k) = sum(RES_TCV(used_els(k,1):used_els(k,2),4));
 end;
 
 Astar = -inv(Lstar)*Rstar;
 Bstar = inv(Lstar); Bstar = Bstar(:,used_coils);
 Bstar = Bstar;

 [Jtcv,Atcv,Rtcv,Ztcv,bp,li,Ic]=tcv_plasma(100,1);
 IP = Jtcv'*Atcv;

 Cstar = 0.5*(Cp+Cm)*IP;

 eArzif = esort(eig(Arzif)); 
 eArzip = esort(eig(Arzip)); 
 eAstar = esort(eig(Astar));

 [eAstar(1:10) eArzip(1:10) eArzif(1:10)] 

 w = logspace(-4,3,30);

 coils = ['OH1';'OH2';'E_1';'E_2';'E_3';'E_4';'E_5';'E_6';'E_7';'E_8';
                      'F_1';'F_2';'F_3';'F_4';'F_5';'F_6';'F_7';'F_8';'G_c'];

 for l=1:2:3
  if(l==1);plot_rz=' R_p ';lr=1;else;plot_rz=' Z_p ';lr=2;end
  for k=1:19
   if(k==1); key ='- = PRO, x = v3(No Ip), + = v3(Ip) '; else ; key = ''; end;

   [Mrzif,Przif] = bode(Arzif,Brzif(:,k),Crzif(lr,:),0,1,w);
   [Mrzip ,Przip ] = bode(Arzip, Brzip( :,k),Crzip( lr,:), 0,1,w);
   [Mstar ,Pstar ] = bode(Astar, Bstar( :,k),Cstar( lr,:), 0,1,w);

   figure(l  );
   plot_MP='Mag: ';
   subplot(5,4,k);cla;
   semilogx(w,20*log10(abs(Mrzip)),'b+',w,20*log10(abs(Mstar)),'r',w,20*log10(abs(Mrzif)),'kx');
   title([key,plot_MP,plot_rz,' vs ',coils(k,1:3)])

   figure(l+1);
   plot_MP='Pha: ';
   subplot(5,4,k);cla;
   semilogx(w,unwrap(pi/180*(Przip)),'b+',w,unwrap(pi/180*(Pstar)),'r',w,unwrap(pi/180*(Przif)),'kx');
   title([key,plot_MP,plot_rz,' vs ',coils(k,1:3)])
  end;
 end;

return

   for k=1:4
    figure(k);
    set(gcf,'PaperUnits','centimeters')
    set(gcf,'PaperPosition', [1, .5, 20.5, 29])
    set(gcf,'PaperUnits','inches')
   end


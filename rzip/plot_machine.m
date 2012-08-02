% To plot the machine described by an RZIP input.
%
%   plot_machine(Machine_name, fig, shot, time);
%
% Inputs
% ======
% Machine_name   the machine description: a local file name
% fig            which figure to plot within
% time           time to get data
% shot           the shot number
%
% if shot ~= 0 then the current filaments are shown. With 10/50/90% boundaries.
%
% J. wainwright Oct 99


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             function   plot_machine(Machine_name, fig, shot, time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALISATION
%%%%%%%%%%%%%%%%

% read in the machine description
 [R_AC, Z_AC, RES_AC, DR_AC, DZ_AC, LA_AC, NT_AC, ...
  R_PS, Z_PS, RES_PS, DR_PS, DZ_PS, LA_PS, ...
  R_BP, Z_BP, TH_BP,                LA_BP, ...
  R_FL, Z_FL,                       LA_FL] = feval(Machine_name);

% set up the figure and grid switch
  figure(fig);  clf;  hold on;

  show_grid_points=1;
  show_passive_names=1;
  show_legend=1;

  col_list=[0.4 0.3 1.0; 0.6 0.3 1.0; 0.8 0.3 1.0; 1.0 0.3 1.0;
            0.6 0.0 0.1; 0.6 0.0 0.3; 0.6 0.0 0.5; 0.6 0.0 0.7];

% Plot the ACTIVE COILS
% ===================== 

% separate active coils into coil groups
% (remember LA_AC names are like 'PF101', 'PF1' = coil name, '01' = coil number in group)

% define a gross list of alternative names (LA_AC(:,1:3))
 [nac,naf,afri,coil_names] = unique_name(LA_AC(:,1:3));

% identify coils in a group
  coil_nums = str2num(LA_AC(:,4:5));

% however a 'coil' could be a compound of coils as well as filaments.
 for k = 1:nac
  coil_range = afri(k,1):afri(k,2);
  coils_in_coil(k,:) = max(coil_nums(coil_range));
 end;

% loops over structures and plot them
 for k = 1:nac
  coil_range = afri(k,1):afri(k,2);
  plot_box(R_AC(coil_range),Z_AC(coil_range),DR_AC(coil_range)/2,DZ_AC(coil_range)/2,'b')
 end;

% display the coils with negative current paths in green
 neg_coils = find(NT_AC<0);
 plot_box(R_AC(neg_coils),Z_AC(neg_coils),DR_AC(neg_coils)/2,DZ_AC(neg_coils)/2,'g')

% loop over coils and label them
 for k = 1:nac
  coil_range=afri(k,1):afri(k,2);
  for l = 1:coils_in_coil(k)
   coil_sub=find(coil_nums(coil_range)==l);
   R_AV = mean(R_AC(coil_range(coil_sub)));
   Z_AV = mean(Z_AC(coil_range(coil_sub)));
   text(R_AV,Z_AV,coil_names(k,1:3),'FontSize',7,'HorizontalAlignment','c','color',[1 1 1],'VerticalAlignment','c')
  end;
 end;


% Plot the PASSIVE STRUCTURES
% ===========================

% passive structure (maybe more than one type)
  [nps, npf, pfri, pass_names] = unique_name(LA_PS(:,1:2));
 
% loop over passive structures and plot them
  for k = 1:nps
   coil_range = pfri(k,1):pfri(k,2);
   plot_box2(R_PS(coil_range),Z_PS(coil_range),DR_PS(coil_range)/2,DZ_PS(coil_range)/2,col_list(k,:),'r')
  end;

% If so desired plot the labels of the passive structures

  if(show_passive_names)

% loop over passive structures and label them by a line at the side
   RMAX = max([R_AC;R_PS]); % maximum R plotted so far
   RMIN = min([R_AC;R_PS]); % minimum R plotted so far
   ZMAX = max([Z_AC;Z_PS]); % maximum Z plotted so far
   ZMIN = min([Z_AC;Z_PS]); % minimum Z plotted so far

   ZBAR = mean([Z_AC;Z_PS]);% average of plotted Z
   RMAX = RMAX*1.1;         % define a position to the right of the plotted area
     DZ = (ZMAX-ZMIN)/20;   % define a size to separate the labels.

   for k = 1:nps;
    coil_range = pfri(k,1):pfri(k,2);
%  find the nearest structure element to the proposed label position
    [i,j] = min(sqrt((R_PS(coil_range)-RMAX).^2+(Z_PS(coil_range)-(ZBAR+(k-1)*DZ)).^2));
    plot([RMAX*0.99 R_PS(coil_range(j))],[ZBAR+(k-1)*DZ Z_PS(coil_range(j))],'Color',get(gca,'XColor'));
    text(RMAX,(ZBAR+(k-1)*DZ),pass_names(k,1:2),'FontSize',8,'HorizontalAlignment','l','Color',col_list(k,:))
   end;
  end;

% Plot the DIAGNOSTICS
% ====================

% flux loops
  for k = 1:length(R_FL)
   plot(R_FL,Z_FL,'bx')
  end;

% bpol probes
  for k = 1:length(R_BP)
   plot(R_BP,Z_BP,'bo')
  end;

% add the normals of the B_POL coils 
  DX = DZ/6;
  for k = 1:length(R_BP);
   RBPN = [R_BP(k)-DX*cos(TH_BP(k)) R_BP(k)+DX*cos(TH_BP(k))];
   ZBPN = [Z_BP(k)-DX*sin(TH_BP(k)) Z_BP(k)+DX*sin(TH_BP(k))];
   plot(RBPN,ZBPN,'b','LineWidth',[3]);
  end;


% If so desired add a legend
% ==========================

  if(show_legend)
   DZ = DZ/5; Fons = 8; 

   plot_box(RMAX,ZMIN,DZ/2,DZ/2,'b');
       atext = text(RMAX+2*DZ,ZMIN,'Active Structure','FontSize',Fons);
   text_size = get(atext,'Extent'); text_sep=1.1*(text_size(4));

   plot_box(RMAX,ZMIN-text_sep,DZ/2,DZ/2,'g');
       atext = text(RMAX+2*DZ,ZMIN-text_sep,'Active but negative path','FontSize',Fons);
   text_size = get(atext,'Extent'); max_extent=text_size(3);

   plot_box(RMAX,ZMIN-2*text_sep,DZ/2,DZ/2,'r');
    text(RMAX+2*DZ,ZMIN-2*text_sep,'Passive Structure','FontSize',Fons);

   plot(RMAX,ZMIN-3*text_sep,'bx');
    text(RMAX+2*DZ,ZMIN-3*text_sep,'Flux loops','FontSize',Fons);

   plot(RMAX,ZMIN-4*text_sep,'bo'); plot([RMAX-DZ RMAX+DZ],[ZMIN-4*text_sep ZMIN-4*text_sep],'b-');
    text(RMAX+2*DZ,ZMIN-4*text_sep,'Bpol Probes','FontSize',Fons);
  end;


% Set up the axis requirements
% turn on zoom
  zoom on
% select appropriate area
 axis([RMIN-2*DZ RMAX+2*DZ+max_extent ZMIN-5*text_sep ZMAX+5*text_sep]);
 axis('equal');
% set up a title
 shot_string='No plasma requested';

% Plot the PLASMA (if required)
% =============================

 if(shot);

  [J_PL,A_PL,R_PL,Z_PL,BETA_P,LI,I_AC] = feval([Machine_name,'_plasma'],shot,time);
   minJ =         0; maxJ = max(J_PL); dJ = (maxJ-minJ)/10;

% plot top 10% in red, 10-50 in mauve 50-90 green 90-100 in yellow
   ind_1 = find(J_PL>9*dJ);
   ind_2 = find(J_PL>5*dJ & J_PL<9*dJ);
   ind_3 = find(J_PL>1*dJ & J_PL<5*dJ); 
   ind_4 = find(J_PL>0    & J_PL<1*dJ); 

   col_PL=[1.0 0.0 0.0 ; 
           1.0 0.6 0.0 ; 
           1.0 0.8 0.0 ; 
           1.0 1.0 0.0 ];

% In absence of DR and DZ for the plasma need to work out best DR and DZ from the long grid.
% The DR and DZ used are the most common values of delta(DR) and delta(DZ) which aren't zero.
   dfr = sort(diff(R_PL)); dfr = dfr(find(dfr));
   dfz = sort(diff(Z_PL)); dfz = dfz(find(dfz));
   [ndis,nnums,rng,nos] = unique_number(dfr); [i,j]=max(rng(:,2)-rng(:,1));   DR_PL=nos(j);
   [ndis,nnums,rng,nos] = unique_number(dfz); [i,j]=max(rng(:,2)-rng(:,1));   DZ_PL=nos(j);

% plot the plasma blocks.
   plot_box2(R_PL(ind_1),Z_PL(ind_1),DR_PL*ones(size(R_PL(ind_1))),DZ_PL*ones(size(R_PL(ind_1))),col_PL(1,:),col_PL(1,:)); hold on
   plot_box2(R_PL(ind_2),Z_PL(ind_2),DR_PL*ones(size(R_PL(ind_2))),DZ_PL*ones(size(R_PL(ind_2))),col_PL(2,:),col_PL(2,:));
   plot_box2(R_PL(ind_3),Z_PL(ind_3),DR_PL*ones(size(R_PL(ind_3))),DZ_PL*ones(size(R_PL(ind_3))),col_PL(3,:),col_PL(3,:));
   plot_box2(R_PL(ind_4),Z_PL(ind_4),DR_PL*ones(size(R_PL(ind_4))),DZ_PL*ones(size(R_PL(ind_4))),col_PL(4,:),col_PL(4,:));

% change title to reflect plasma presence.
   shot_string=['Plasma Shot Number: ',num2str(shot,5),' at time ',num2str(time),'s'];
 
 end;

% display a title.

  eval(['title(''Machine ',Machine_name,' : ',shot_string,''')']);

 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES required by Plot_machine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%                                         plot_box2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a box of dimensions dr, dz at r,z with face colour colf and border colour colb.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function plot_box2(r,z,dr,dz,colb,colf);

 nr=length(r);
 for k=1:nr 
  x=[r(k)-dr(k) r(k)-dr(k) r(k)+dr(k) r(k)+dr(k) r(k)-dr(k)];
  y=[z(k)+dz(k) z(k)-dz(k) z(k)-dz(k) z(k)+dz(k) z(k)+dz(k)];
  fill(x,y,colf);
  x=[r(k)-dr(k) r(k)-dr(k) r(k)+dr(k) r(k)+dr(k) r(k)-dr(k)];
  y=[z(k)+dz(k) z(k)-dz(k) z(k)-dz(k) z(k)+dz(k) z(k)+dz(k)];
  plot(x,y,'color',colb);
 end;

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
%                                       unique_number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a file which give a column vector of numbers, counts the different numbers and defines the
% range in the list of each name (index range) and total length of vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [nn,nf,index_range,names]=unique_number(numbers);

% define total length and create a row of the names
  nf=size(numbers,1); 

% a potential error is if a name is made of all the same characters!
% To avoid this scenario add a dummy separation variable
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


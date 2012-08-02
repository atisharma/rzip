% This file takes in some basic passive data and creates and input for tcv_tiles_plasma.m
% we interpolate np points extra in the structure as we are replacing with
% rectangles and its too chunky otherwise 

  function make_tiles(tile_res,tile_range);

% load in the old data
  load tiles_in
  tile_width=.01;
%  if(tile_res==999);tile_res=13.75e-6;end;

  t_r_t=[]; t_z_t=[];

% stupid bloody definition of tiles
%  interpolate np_array points between certain lengths 
   np_array=[6 6 6 3 6 6 3 6 6 3];
   
   for k=1:10
    dr_t=(t_r(k+1)-t_r(k))/np_array(k);
    dz_t=(t_z(k+1)-t_z(k))/np_array(k);
    dummy_r=ones(1,np_array(k))*dr_t; dummy_r(1)=0;
    dummy_z=ones(1,np_array(k))*dz_t; dummy_z(1)=0;
    t_r_t=[t_r_t [t_r(k)+cumsum(dummy_r)]];
    t_z_t=[t_z_t [t_z(k)+cumsum(dummy_z)]];
   end;
 
   t_r_t=[t_r_t t_r(11:4:38)];
   t_z_t=[t_z_t t_z(11:4:38)];

   np_array=2;
   dr_t=(t_r(41)-t_r(40))/np_array;
   dz_t=(t_z(41)-t_z(40))/np_array;
   dummy_r=ones(1,np_array)*dr_t; dummy_r(1)=0;
   dummy_z=ones(1,np_array)*dz_t; dummy_z(1)=0;
   t_r_t=[t_r_t [t_r(40)+cumsum(dummy_r)]];
   t_z_t=[t_z_t [t_z(40)+cumsum(dummy_z)]];


% will increase the number of points by interpolation of np points

  [RR,ZZ,XX,YY,NN]=add_points_closed(t_r_t',t_z_t',0);
  arc_length=sqrt(diff(RR).^2+diff(ZZ).^2);
  for k=1:NN;
   DR=RR(k+1)-RR(k);
   DZ=ZZ(k+1)-ZZ(k);
   if(DR~=0)
    grad = DZ/DR;
   else
    grad=1e20;
   end;
   if(abs(grad)>1);
    DR_TL(k)=tile_width;
    DZ_TL(k)=arc_length(k);
   else;
    DR_TL(k)=arc_length(k);
    DZ_TL(k)=tile_width;
   end;
  end;

  R_TL=XX;
  Z_TL=YY;
  RES_TL=2*pi*tile_res*R_TL./(DR_TL.*DZ_TL);

  save tiles.mat R_TL Z_TL RES_TL DR_TL DZ_TL tile_range 

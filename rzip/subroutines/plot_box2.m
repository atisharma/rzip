function plot_box(r,z,dr,dz,colb,colf);

 nr=length(r);
% hold on;
 for k=1:nr 
  x=[r(k)-dr(k) r(k)-dr(k) r(k)+dr(k) r(k)+dr(k) r(k)-dr(k)];
  y=[z(k)+dz(k) z(k)-dz(k) z(k)-dz(k) z(k)+dz(k) z(k)+dz(k)];
  fill(x,y,colf);
  x=[r(k)-dr(k) r(k)-dr(k) r(k)+dr(k) r(k)+dr(k) r(k)-dr(k)];
  y=[z(k)+dz(k) z(k)-dz(k) z(k)-dz(k) z(k)+dz(k) z(k)+dz(k)];
  plot(x,y,colb);
 end;

return;

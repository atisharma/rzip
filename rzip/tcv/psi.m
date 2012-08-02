
% Get current profile
[J,A,R,Z,bp,li,I_AC]=tcv_plasma(100,1);

load /home/ps/sunfs1/d4/ati/matlab/model/psi13333
load tcv_plasma_full_v2
psi_ext= [I_AC.'*mut_struct_grid(1:19,:)].';
psi_self2 = PSI(:) - psi_ext;

% Define plasma radius,
% and other plasma properties/parameters.
nz=find(J ~=0);
Ip = J.'*A;
R0 = sum(J.*A.*R)/Ip;
Z0 = sum(J.*A.*Z)/Ip;

dR = R-R0;
dZ = Z-Z0;
lambda = bp + li/2 -1;
mu0 = 4*pi*1e-7;
a = (max(R(nz))-min(R(nz)))/2;
b = (max(Z(nz))-min(Z(nz)))/2;
kappa = b/a;
rho = sqrt(dR.*dR + (dZ.*dZ)/kappa/kappa);
omega = atan2(dZ/kappa,dR);
gamma = (log(8*R0/a/sqrt(kappa)) + bp + li/2 -3/2);

% Calculate PSI, plasma self-flux, and derivatives

psi_self = (R0*Ip*mu0).*(log(8*R0./rho) - 2) + mu0*Ip/2*(log(8*R0./rho) - 1 + a*a./rho./rho*(lambda+0.5)).*rho.*cos(omega);
dpsi_dR = mu0 * R0 * Ip .*(log(8.*R0./rho)-1) + mu0/2*Ip./R0;

% Plasma self-inductance
Lp = 1/Ip/Ip*sum(J(nz).*psi_self2(nz).*A(nz));
Lp_old = mu0*R0*(log(8*R0/a/sqrt(kappa)) - 2 + li/2);
Lp_circ = mu0*R0*(log(8*R0/a) - 2 + li/2);
dLp_dR = 1/Ip/Ip*sum(J(nz).*dpsi_dR(nz).*A(nz));
M34 = mu0*gamma;
M43 = mu0*(log(8*R0/a/sqrt(kappa)) - 1 + li/2);

%------------------------------------------------------

% plot profiles
%figure(1)
%clf
dd=[0:10000:max(J)];                                    
%plot3(dR(nz),dZ(nz),J(nz),'ro')
hold on 
%plot3(zeros(size(dd)),zeros(size(dd)),dd,'g.')
%plot3(dR,Z,zeros(size(J)),'k.')
axis('equal')
title(' Plasma current profile')
xlabel('R - R_0')
ylabel('Z - Z_0')
zlabel('J')

figure(2)
clf
%plot3(dR,dZ,psi_self,'b.')
%grid on
%hold on
%plot3(dR(nz),dZ(nz),psi_self(nz),'r.')
%title('Psi; plasma self-flux')
%xlabel('R - R_0')
%ylabel('Z - Z_0')
%zlabel('Psi')



Ip
kappa
Lp
Lp_old
ratio = Lp/Lp_old

dLp_dR
M34
M43



A = dR(1:28);
B = dZ(1:28:end);
Cp2 = reshape(psi_self2,length(A),length(B));
Cp = reshape(psi_self,length(A),length(B));
Cj = reshape(J,length(A),length(B));
Cpe= reshape(psi_ext,length(A),length(B));

figure(3)
clf
contour(A,B,Cp2',10,'r')
hold on
%contour(A,B,Cj',10,'b')
%contour(A,B,PSI',10,'k')
%contour(A,B,Cpe',10,'c')
contour(A,B,Cp',10,'b')
grid on
axis equal
%axis([min(dR(nz)) max(dR(nz)) min(dZ(nz)) max(dZ(nz))])
xlabel('dR')
ylabel('dZ')
title('J and Psi contours')



figure(2)
mesh(A,B,Cp2')
hold on
%mesh(A,B,Cp')
mesh(A,B,Cpe')
mesh(A,B,PSI')
xlabel('dR')
ylabel('dZ')

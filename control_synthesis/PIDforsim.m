
% Save discrete time controller
%
% expects [Bpol; Psi; Icoil(1:19)]


PID = ss(Ak,Bk,-Ck,-Dk,1e-4);
PID = PID(1:18,:);
%Gps_dt = c2d(Gps,1e-4,'tustin');

%display('Including power supply');
%PID = Gps_dt*PID;

AK = PID.A;
BK = PID.B;
CK = PID.C;
DK = PID.D;

save PID_controller.mat AK BK CK DK 

%clear all
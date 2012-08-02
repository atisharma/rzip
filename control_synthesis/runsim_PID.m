% Load controller, model, and matrices needed by simulation.
% compares to PROTEUS simulation data

%load simmat
%load rzipsim y

% assume Hinf simulation already run
y_Hinf = y;clear y
y_Hinf(:,5) = y_Hinf(:,5)/1e4;

load PID_controller

delta_t = 1e-4; % fixed

%open_system('conttest_discrete_PID')
sim('conttest_discrete_PID')

ylabels = ['p_{vert}    ';'tri_{out}   ';'tri_{in}    ';'zI_{p}      ';'I_{p}/10^{4}';'\psi_{R}    '];

y(:,5)=y(:,5)/1e4;
ref(:,5)=ref(:,5)/1e4;

figure(1)
clf
for n=1:5
   m=n; if m==1
         m=4;
      elseif m==4
         m=1;
      end
   subplot(5,1,m)
%   plot(T,ref(:,n),'r',T,y(:,n),'b',t,destim(:,n),'g')
%plot(T,ref(:,n),'r',T,y(:,n),'g')
plot(T,ref(:,n),'r')
hold on
plot(T(1:length(y)),y(:,n), T(1:length(y_Hinf)),y_Hinf(:,n))
   ylabel(ylabels(n,:))
   if m==5
      xlabel('T [s]')
   end
   if m==1
      title('Tracking of control parameters with standard PID controller, H_{\infty} controller')
   end
end
%orient landscape; print
%print -depsc2 D:/latex/thesis/figures/tracking_PID_alt.eps


figure(3)
clf
plot(T,u)
ylabel('P/S demand [V]')
xlabel('T [s]')
title('PF coil P/S demand signal from standard PID controller')
%orient landscape; print 
%print -depsc2 D:/latex/thesis/figures/udemand_PID_alt.eps

y_PID = y;
u_PID = u;

save rzipsim_PID T ref y* u_PID
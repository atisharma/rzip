% Load controller, model, and matrices needed by simulation.

%open_system('conttest_discrete')
sim('conttest_discrete')

ylabels = ['p_{vert} ';'tri_{out}';'tri_{in} ';'zI_{p}   ';'I_{p}    ';'\psi_{R} '];

figure(1)
clf
for n=1:5
   m=n; if m==1
         m=4;
      elseif m==4
         m=1;
      end
   subplot(5,1,m)
   plot(T,ref(:,n),'r',T,y(:,n),'b')
%   plot(T,ref(:,n),'r',T,y_PID(:,n),'b', T,y(:,n),'g')
   ylabel(ylabels(n,:))
   if m==5
      xlabel('T [s]')
   end
   if m==1
      title('Tracking of control parameters with H_{\infty} controllers')
   end
end


figure(2)
clf
plot(T,u)
ylabel('P/S demand [V]')
xlabel('T [s]')
title('PF coil P/S demand signal from H_{\infty} controller')

%save rzipsim T ref y* u 
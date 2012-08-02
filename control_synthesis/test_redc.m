% A Sharma 2001

clear all
gamma = inf;
load names
% load tokamak model
load rmod
% load power supply model & limits
load ps

% Include p/s
G=G*Gps;

[Wo,Wi,Gw] = weights(G);
%Gw = Wo*G*Wi;

% l = [order;gamma;eig]
l=[];

for order = (5:35)
   disp(['Finding controller of order ' num2str(order) '...'])
   [Kwred, S, bound] = redc(Gw, order);
   Kred = Wi*Kwred*Wo;
   Hwred = closeloop2(Gw,Kwred);
   Hred = closeloop2(G,Kred);
   er = max(real(esort(eig(Hwred.A))));
   n = norm(Hwred,inf);
   l = [l [order;n;bound;er]];
   disp(['inf norm:' num2str(n) ' bound:' num2str(bound) ', max eig: ' num2str(er)] )
end

%save reduced_list l

figure(1);clf
s=diag(S);
semilogy(1:length(s),s);
title('sigma')

figure(2);clf
s = find(l(4,:)<0);
plot(l(1,s),l(2,s),'b*',l(1,s),l(3,s),'r+')
grid on
xlabel('{order of Kred}')
ylabel('{replaceme}')
legend('    {replaceme}        ','upper bound')
title('||\cdot||_{\infty} of closed loop (actual and bound) against order of reduced controller ')
%print -depsc2 D:/latex/IEEEpaper/figures/controller_order_alt.eps

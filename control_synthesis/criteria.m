% Make singular value plots to check 
% loop shapes in various ways
% Plot all relevant figures
w = logspace(-2,6,50);


figure(3)
clf
grid on
sigma(Gw, w);
title(['Singular value plot of weighted plant (including P/S)'])
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')
if exist('gamma')
   text(1e4,150,['\gamma = ' num2str(gamma)])
end
%orient landscape; print



% Loop shaping relevant figures
figure(5)
clf
grid on
sigma((G)*(Wi*Kw*Wo), w)
title(' \sigma_{min}(GK) >> 1 for low \omega, \sigma_{max}(GK) << 1 for high \omega')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')
%orient landscape; print

figure(6)
clf
grid on
sigma((Wi*Kw*Wo)*(G), w)
title(' \sigma_{min}(KG) >> 1 for low \omega, \sigma_{max}(KG) << 1 for high \omega')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')
%orient landscape; print

figure(7)
grid on
sigma(Wi*Kw*Wo, w)
title(' \sigma_{min}(K) >> 1 for low \omega, \sigma_{max}(K) <= M for high \omega')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')
%orient landscape; print


% Produce figure for controller reduction
figure(4)
clf
hkv=hksv(-K.A,-K.B,K.C);
plot(hkv,'+')
title('Hankel singular value plot of controller K')
ylabel('HSV')

figure(1)
clf
[sig2,w] = sigma(Hw,w);
semilogx(w,sig2,'b')
grid on
title(['Singular value plot of closed loop (weighted) plant, (\gamma_{opt} = ' num2str((gamma)) ')'])
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value')


figure(2)
clf
grid on
sigma(G,w);
title(['Singular value plot of unweighted plant (including P/S)'])
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')

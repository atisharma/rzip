% Make singular value plots to check 
% loop shapes in various ways
% Plot all relevant figures
w = logspace(-2,6,50);

% Loop shaping relevant figures
figure(5)
clf
grid on
sigma((G)*(Kred), w)
hold on
sigma((G)*(K), w)
title(' \sigma_{min}(G_{W}K) >> 1 for low \omega, \sigma_{max}(G_{W}K) << 1 for high \omega')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')
%orient landscape; print

figure(6)
clf
grid on
sigma((Kred)*(G), w)
hold on
sigma((K)*(G), w)
title(' \sigma_{min}(KG_{W}) >> 1 for low \omega, \sigma_{max}(KG_{W}) << 1 for high \omega')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')
%orient landscape; print

figure(7)
clf
grid on
sigma(Kred, w)
hold on
sigma(K,w)
title(' \sigma_{min}(K) >> 1 for low \omega, \sigma_{max}(K) <= M for high \omega')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')
%orient landscape; print


figure(1)
clf
[sig2,w] = sigma(Hwred,w);
[sigopt,w] = sigma(Hw,w);
semilogx(w,sig2,'r--')
hold on
semilogx(w,sigopt,'b')
grid on
title(['Singular value plot of closed loop (weighted) plant, (\epsilon = 'num2str((1/norm(Hwred,inf))) ')'])
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value')



figure(8)
clf
grid on
sigma(Hwred,w,'r--');
hold on
sigma(Hw,w,'b');
title(['Singular value plot of closed loop (weighted) plant, (\epsilon = 'num2str((1/norm(Hwred,inf))) ')'])
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')


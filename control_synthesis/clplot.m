load testhinfcont


w = logspace(-1,6,300);
svH  = sigma(Hw, w);
svHr = sigma(Hwred, w);

figure(1);clf
semilogx(w, svH, 'b');
grid on
hold on
semilogx(w, svHr, 'r--')

title('Singular value plot of closed loop (weighted) plant')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value')


figure(2);clf
grid on
svK = sigma(K,w);
svKr = sigma(Kred,w);
semilogx(w,20*log10(svK),'b');
hold on
semilogx(w,20*log10(svKr),'r--');
axis([1e-1 1e6 -80, 120])
grid on

title('Singular values of the controllers')
xlabel('Frequency (rad s^{-1})')
ylabel('Singular value (dB)')


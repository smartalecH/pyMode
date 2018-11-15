
clf

subplot(3,1,1)
semilogy(lambdas/1e-6, dneffte, '+-');
hold on
semilogy(lambdas/1e-6, dnefftm, 'r+-');
title 'Relative error in effective indices'
legend('TE','TM','Location','NorthWest')
xlim([lambdas(1) lambdas(end)] / 1e-6)

subplot(3,1,2)
semilogy(lambdas/1e-6, dTEez, '+-', 'LineWidth', 2);
hold on
semilogy(lambdas/1e-6, dTMer, 'r+-', 'LineWidth', 2);
semilogy(lambdas/1e-6, dTEer, 'o-');
semilogy(lambdas/1e-6, dTEep, 'o:');
semilogy(lambdas/1e-6, dTMez, 'ro--');
semilogy(lambdas/1e-6, dTMep, 'ro:');
legend('TE E_z','TM E_\rho','Location','NorthWest')
title 'Error measure for electric-field components'
xlim([lambdas(1) lambdas(end)] / 1e-6)
    
subplot(3,1,3)
semilogy(lambdas/1e-6, dTEhr, '+-', 'LineWidth', 2);
hold on
semilogy(lambdas/1e-6, dTMhz, 'r+-', 'LineWidth', 2);
semilogy(lambdas/1e-6, dTEhz, 'o--');
semilogy(lambdas/1e-6, dTEhp, 'o:');
semilogy(lambdas/1e-6, dTMhr, 'ro-');
semilogy(lambdas/1e-6, dTMhp, 'ro:');
legend('TE H_\rho','TM H_z','Location','NorthWest')
title 'Error measure for magnetic-field components'
xlim([lambdas(1) lambdas(end)] / 1e-6)

xlabel 'Wavelength / um'

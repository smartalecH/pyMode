figure(1); clf

%load doit01.mat
%semilogy(PARs, 1e4*Alphas', 'x-')
%hold on
%load doit02.mat
semilogy(PARs, 1e-2*Alphas', 'x--')

ylabel 'Alpha [dB/cm]'

xlim([1.0 1.6]*1e-6)
ylim([1e-3 1e3])

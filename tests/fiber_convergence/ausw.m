
figure(1); clf

set(gcf, 'PaperType', 'A4')
set(gcf, 'PaperOrientation', 'portrait')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [1 1 12 7.5])
set(gca, 'FontSize', 10)

load case13_1_dnonzero
relerrs = abs(neffs(1,:) - neffreal(1)) / neffreal(1);
loglog(XMAX./(Ns-1), relerrs, '^-', 'LineWidth', .75, 'MarkerSize', 5) % TE
relerrs = abs(neffs(2,:) - neffreal(2)) / neffreal(2);
hold on
loglog(XMAX./(Ns-1), relerrs, 'o-', 'LineWidth', .75,  'MarkerSize', 5) % HE21

load case13_2_dnonzero
relerrs = abs(neffs(1,:) - neffreal(1)) / neffreal(1);
loglog(XMAX./(Ns-1), relerrs, 'x-', 'LineWidth', .75,  'MarkerSize', 5) % HE11 = fund
relerrs = abs(neffs(2,:) - neffreal(2)) / neffreal(2);
loglog(XMAX./(Ns-1), relerrs, 's-', 'LineWidth', .75,  'MarkerSize', 5) % EH11

load case13_3_dnonzero
relerrs = abs(neffs(1,:) - neffreal(1)) / neffreal(1);
loglog(XMAX./(Ns-1), relerrs, 'x-') % HE11 = fund
                                    % identical curve!
relerrs = abs(neffs(2,:) - neffreal(2)) / neffreal(2);
loglog(XMAX./(Ns-1), relerrs, 'kx-') % EH11
                                     % identical curve!

load case13_4_dnonzero
relerrs = abs(neffs(1,:) - neffreal(1)) / neffreal(1);
loglog(XMAX./(Ns-1), relerrs, 'o-', 'LineWidth', .75,  'MarkerSize', 5) % HE21
relerrs = abs(neffs(2,:) - neffreal(2)) / neffreal(2);
loglog(XMAX./(Ns-1), relerrs, 'v-', 'LineWidth', .75,  'MarkerSize', 5) % TM


%load case13_2_dzero
%relerrs = abs(neffs(1,:) - neffreal(1)) / neffreal(1);
%loglog(XMAX./(Ns-1), relerrs, 'x--', 'LineWidth', .75, 'MarkerSize', 5)


%load case13_3_dzero
%relerrs = abs(neffs(1,:) - neffreal(1)) / neffreal(1);
%loglog(XMAX./(Ns-1), relerrs, 'x-', 'LineWidth', 2)
% identical curve!


xlabel 'Grid-point distance [um]'
ylabel 'Relative error in n_{eff}'
xlim([1e-2 4e-1])

print -deps2 convergence.eps

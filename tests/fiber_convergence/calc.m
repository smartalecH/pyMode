
Ns = logspace(log10(40),log10(1000),20);
%Ns = linspace(50,700,25);
Ns = floor(Ns);

ms = [ ms ' ' symm ];
ls = { 'x-', 'ro-' };

neffs = nan * ones(length(neffreal), length(Ns));

XMAX = 12.99101;

for Nc = 1 : length(Ns)
    com = sprintf('%s -U 0:%d:%e -V 0:%d:%e', ms, Ns(Nc), XMAX, Ns(Nc), XMAX);
    m = wgms3d_run(com);
    for i = 1 : length(neffreal)
        neffs(i,Nc) = m{i}.neff;
    end
  
    cfig(4); clf

    for i = 1 : length(neffreal)
        relerrs = abs(neffs(i,:) - neffreal(i)) / neffreal(i);
        loglog(XMAX./(Ns-1), relerrs, ls{i})
        hold on
    end
    xlabel 'Grid-point distance [um]'
    ylabel 'Relative error in n_{eff}'
    %    xlim(XMAX./([max(Ns) min(Ns)]-1))
    xlim([1e-2 XMAX/(min(Ns)-1)])
    drawnow
end

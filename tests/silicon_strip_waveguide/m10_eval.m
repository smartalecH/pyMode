
clear cf
clear nes ngs Ds

ls = linspace(lambdas(1), lambdas(end), 101);
c = 299792458;

cfig(10); clf
set(gcf, 'PaperType', 'A4')
set(gcf, 'PaperOrientation', 'portrait')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [ .5 .5 get(gcf,'PaperSize')-1 ]);
set(gca, 'FontSize', 10)


titles = { 'quasi-TE', 'quasi-TM' };

for i = 1 : 2
    subplot(4,2,1+(i-1))
    plot(lambdas/1e-6, neffs(i,:), 'x')
    title(titles{i})
    ft_ = fittype('poly5');
    fo_ = fitoptions('method','LinearLeastSquares','Normalize','on');
    %fo_ = fitoptions('method','CubicSplineInterpolant','Normalize','on');
    cf{i} = fit(lambdas', neffs(i,:)', ft_, fo_);
    
    for l = 1 : length(ls)
        nes(i,l) = feval(cf{i}, ls(l));
        % Das Buch, (4.139)
        [nel, nell] = differentiate(cf{i}, ls(l));
        ngs(i,l) = nes(i,l) - ls(l)*nel;
        % Das Buch, (4.143)
        D(i,l)   = -(ls(l)/c)*nell;
    end
    
    hold on
    plot(ls/1e-6, nes(i,:))
    xlim([min(ls) max(ls)]/1e-6)
    ylabel 'neff'
    
    subplot(4,2,3+(i-1))
    plot(ls/1e-6, ngs(i,:))
    xlim([min(ls) max(ls)]/1e-6)
    ylabel 'ng'
    
    subplot(4,2,5+(i-1))
    % 1 ps/(nm*km) = 1e-12 / 1e-9 / 1e3 s/m^2 = 1e-6 s/m^2
    % 1 fs/(nm*cm) = 1e-15 / 1e-9 / 1e-2 s/m^2 = 1e-4 s/m^2
    plot(ls/1e-6, D(i,:)/1e-4)
    ylabel 'D [fs/(nm*cm)]'
    xlim([min(ls) max(ls)]/1e-6)
    
    % Das Buch, (4.146)
    subplot(4,2,7+(i-1))
    ft_ = fittype('poly5');
    fo_ = fitoptions('method','LinearLeastSquares','Normalize','on');
    %fo_ = fitoptions('method','CubicSplineInterpolant','Normalize','on');
    Dfit = fit(ls', D(i,:)', ft_, fo_);
    % 1 ps/(nm^2*km) = 1e-12 / 1e-18 / 1e3 s/m^2 = 1e3 s/m^2
    % 1 fs/(nm^2*cm) = 1e-15 / 1e-18 / 1e-2 s/m^2 = 1e5 s/m^2
    plot(ls/1e-6, differentiate(Dfit, ls)/1e5)
    ylabel 'S [fs/(nm^2*cm)]'
    xlim([min(ls) max(ls)]/1e-6)
    xlabel 'Wavelength [nm]'
    
end

print -dps2c m10_eval.ps

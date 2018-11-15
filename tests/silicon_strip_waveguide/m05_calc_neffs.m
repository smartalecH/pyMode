
% called from main.m

% set W, H, RC, lambdas

xx = wgms3d_grid_make([   -1.2e-6  -W/2-100e-9  -W/2-25e-9  -W/2+25e-9  -W/2+100e-9     0], ...
                      [   6.01e-9      6.01e-9     2.01e-9     2.01e-9      6.01e-9  6.01e-9]);
xx = xx - xx(end);

yy = wgms3d_grid_make([   -2.0e-6  -H/2-100e-9  -H/2-25e-9  -H/2+25e-9  -H/2+100e-9     0], ...
                      [   6.01e-9      6.01e-9     2.01e-9     2.01e-9      6.01e-9  6.01e-9]);
yy = yy - yy(end);
yy = [ yy(1:end-1) -fliplr(yy) ];
yy = yy + H/2;

ms = 'LD_LIBRARY_PATH=. ~/wgms3d/current/wgms3d';

wgms3d_grid_write('xx.txt', xx);
wgms3d_grid_write('yy.txt', yy);

clear modes
neffs = nan(2, length(lambdas));

for lambdac = 1 : length(lambdas)
    lambda = lambdas(lambdac);
    fn = 'mgp.mgp';
    wgms3d_mgp_rib_waveguide(fn, W, H, H, RC, 1.0, get_nSi(lambda/1e-6), get_nSiO2(lambda/1e-6), []);

    % quasi-TE mode:
    modes = wgms3d_run(sprintf(...
        ['%s -g %s -l %e -U xx.txt -V yy.txt -n 1'], ...
        ms, fn, lambda));
    neffs(1,lambdac) = modes{1}.neff;
    cfig(10); clf
    wgms3d_plot_ht(0, 'Geometry', 'mgp.mgp', 'LogContours', 3)
    drawnow
        
    % quasi-TM mode:
    modes = wgms3d_run(sprintf(...
        ['%s -g %s -l %e -U xx.txt -V yy.txt -n 1 -M e'], ...
        ms, fn, lambda));
    neffs(2,lambdac) = modes{1}.neff;
    cfig(11); clf
    wgms3d_plot_ht(0, 'Geometry', 'mgp.mgp', 'LogContours', 3)
    drawnow
end

wgms3d_clean


function com = nguyen_prepfunc (W, nsearch)
    
    lambda0 = 1.55e-6;
    
    H = 0.205e-6;
    D = 0.015e-6;
    RC = 0.007e-6;
    wgms3d_mgp_rib_waveguide('nguyen.mgp', W, H, D, RC, 1.0, 3.4797, 1.444);

    xx = wgms3d_grid_make([     0e-6  W/2-100e-9   W/2-25e-9   W/2+25e-9   W/2+100e-9   5.0e-6], ...
                          [  8.01e-9     8.01e-9     0.51e-9     0.51e-9      8.01e-9  8.01e-9]);
    
    yy = wgms3d_grid_make([  -1.5e-6   H-100e-9     H-25e-9   H+25e-9  H+100e-9   1.3e-6], ...
                          [  8.01e-9    8.01e-9     0.51e-9   0.51e-9   8.01e-9  8.01e-9]);

    wgms3d_grid_write('xx.txt', xx);
    wgms3d_grid_write('yy.txt', yy);

    ms = 'LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -U xx.txt -V yy.txt -M w -P e:30:0.1';
    
    com = sprintf('%s -g nguyen.mgp -n 2 -E -F -l %e -s %e', ...
                  ms, lambda0, abs(mean(nsearch)));
    

% To get a good agreement between the approximate and exact results, it is
% important to make the grid fine enough in the upper-cladding region, since
% the field decays rapidly there. Also, the two waveguide corners must be
% resolved well, since there is a near-singularity in the electric
% field. Therefore there must be many grid points for a stable computation
% of the integral in these regions.

W  = 700e-9;
H  = 220e-9;
D  =  70e-9;
RC =   5e-9;

xx = wgms3d_grid_make([     0e-6  W/2-100e-9   W/2-25e-9   W/2+25e-9   W/2+100e-9   1.5e-6], ...
                      [  1.89e-9     1.89e-9     0.57e-9     0.57e-9      8.01e-9  8.01e-9]);
    
yy = wgms3d_grid_make([  -1.5e-6   80e-9     100e-9       H+25e-9  H+100e-9  H+200e-9 H+300e-9 1.3e-6], ...
                      [  8.01e-9    8.01e-9     0.57e-9   0.57e-9   1.89e-9  1.89e-9  8.01e-9  8.01e-9]);
length(xx)
length(yy)

wgms3d_grid_write('xx.txt', xx);
wgms3d_grid_write('yy.txt', yy);

ms = 'LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -U xx.txt -V yy.txt -l 1.55e-6';

alpha = .2303 * 1.0; % 1 dB/m

% Calculate loss-less mode:

wgms3d_mgp_rib_waveguide('lossless.mgp', W, H, D, RC, 1.0, 3.5, 1.5);
com = sprintf('%s -g lossless.mgp -n 1 -e -E -F -G', ms);
m1 = wgms3d_run(com);
m1{1} = wgms3d_load_mode_field(m1{1});
wgms3d_plot_et(m1{1},'Geometry','lossless.mgp','QuiverGrid',0)
axis equal
drawnow

% Calculate loss from lossless mode via perturbation theory [e.g., Snyder /
% Love, Eq. (18-71)]: (assume the top cladding is lossy)

Z0 = 376.7303135;
mask = 1.0 * (m1{1}.epsis == (1.0)^2);
m1{1}.N = wgms3d_modeproduct(m1{1}, m1{1});
aa = abs(m1{1}.er).^2 + abs(m1{1}.ez).^2 + abs(m1{1}.ep).^2;
eta = wgms3d_int(aa .* mask, m1{1}.r, m1{1}.z) / (2 * m1{1}.N * Z0)

disp(sprintf('approximate: %e dB/m', eta*alpha / .2303))

% -----------------------------------------------------------------------

% Calculate lossy mode

% Note that there exist several conflicting conventions regarding the
% relation between the complex refractive index and the "extinction
% coefficient" (or "attenuation index") kappa. We here use the notation of
% Born / Wolf, Chapter 14.1, such that the complex refractive index may be
% decomposed in real and imaginary parts according to
%
%   n_complex = n_real * (1 + i * kappa)
%
% Under this convention, kappa is related to the power-absorption
% coefficient alpha through
%
%   kappa = alpha / (2 * k0 * n_real),
%
% where k0 is the free-space wavenumber.

k0 = 2*pi/1.55e-6;
kappa = alpha / (2 * k0 * 1.0);

wgms3d_mgp_rib_waveguide('lossy.mgp', W, H, D, RC, 1.0*(1+kappa*i), 3.5*(1+0*i), 1.5*(1+0*i));
com = sprintf('%s -g lossy.mgp -n 1 -e', ms);
m2 = wgms3d_run(com);
m2{1} = wgms3d_load_mode_field(m2{1});

disp(sprintf('exact:       %e dB/m', m2{1}.alphadbm))


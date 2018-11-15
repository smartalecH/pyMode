
% This script takes a fixed waveguide geometry and calculates the two
% fundamental modes (quasi-TE and quasi-TM) at several wavelengths. It then
% calculates and plots the error that is made by the semi-vectorial (SV)
% approximation with respect to the exact full-vectorial one.
%
% As expected, the error made using the SV approximation becomes smaller
% when the wavelength decreases, i.e., when the mode is more strongly
% guided.

% Note that the precise form of the computed curves for the E-field error
% depends significantly on the grid chosen, i.e. on how fine you sample the
% region around the dielectric "corner" where the E field has a more or less
% pronounced singularity.

ms = 'LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -E -F -G -H -U -4e-6:300:0 -V -.995e-6:251:2.505e-6';

wgms3d_mgp_rib_waveguide('rib.mgp', 1.5e-6, 1.5e-6, 0.75e-6, 50e-9, ...
                         1.0, 3.5, 1.5);

lambdas = linspace(1.0e-6, 4.5e-6, 50);

dneffte = nan(1, length(lambdas));
dnefftm = dneffte;

dTEer = dneffte;
dTEez = dTEer;
dTEep = dTEer;
dTMer = dTEer;
dTMez = dTEer;
dTMep = dTEer;
dTEhr = dTEer;
dTEhz = dTEer;
dTEhp = dTEer;
dTMhr = dTEer;
dTMhz = dTEer;
dTMhp = dTEer;

figure(1); clf
drawnow

for lambdac = 1 : length(lambdas)
        
    lambda = lambdas(lambdac);

    % Calculate full-vectorial modes
    
    com = sprintf('%s -l %e -g rib.mgp -n 1', ms, lambda);
    m = wgms3d_run(com);
    FVTE = wgms3d_load_mode_field(m{1});
    
    com = sprintf('%s -l %e -g rib.mgp -n 1 -M e', ms, lambda);
    m = wgms3d_run(com);
    FVTM = wgms3d_load_mode_field(m{1});

    % Now calculate semi-vectorial modes

    com = sprintf('%s -l %e -g rib.mgp -n 1 -u', ms, lambda);
    m = wgms3d_run(com);
    SVTE = wgms3d_load_mode_field(m{1});
    
    com = sprintf('%s -l %e -g rib.mgp -n 1 -v -M e', ms, lambda);
    m = wgms3d_run(com);
    SVTM = wgms3d_load_mode_field(m{1});

    % Compare effective indices
    
    dneffte(lambdac) = abs(FVTE.neff - SVTE.neff) / FVTE.neff;
    dnefftm(lambdac) = abs(FVTM.neff - SVTM.neff) / FVTM.neff;
    
    % Compare mode fields
    
    FVTE = wgms3d_normalize_mode_field(FVTE);
    FVTM = wgms3d_normalize_mode_field(FVTM);
    SVTE = wgms3d_normalize_mode_field(SVTE);
    SVTM = wgms3d_normalize_mode_field(SVTM);

    FVTEd = wgms3d_int(abs(FVTE.er).^2+abs(FVTE.ez).^2+abs(FVTE.ep).^2, FVTE.r, FVTE.z);
    FVTMd = wgms3d_int(abs(FVTM.er).^2+abs(FVTM.ez).^2+abs(FVTM.ep).^2, FVTM.r, FVTM.z);
    dTEer(lambdac) = wgms3d_int((abs(FVTE.er)-abs(SVTE.er)).^2, FVTE.r, FVTE.z) / FVTEd;
    dTEez(lambdac) = wgms3d_int((abs(FVTE.ez)-abs(SVTE.ez)).^2, FVTE.r, FVTE.z) / FVTEd;
    dTEep(lambdac) = wgms3d_int((abs(FVTE.ep)-abs(SVTE.ep)).^2, FVTE.r, FVTE.z) / FVTEd;
    dTMer(lambdac) = wgms3d_int((abs(FVTM.er)-abs(SVTM.er)).^2, FVTM.r, FVTM.z) / FVTMd;
    dTMez(lambdac) = wgms3d_int((abs(FVTM.ez)-abs(SVTM.ez)).^2, FVTM.r, FVTM.z) / FVTMd;
    dTMep(lambdac) = wgms3d_int((abs(FVTM.ep)-abs(SVTM.ep)).^2, FVTM.r, FVTM.z) / FVTMd;

    FVTEd = wgms3d_int(abs(FVTE.hr).^2+abs(FVTE.hz).^2+abs(FVTE.hp).^2, FVTE.r, FVTE.z);
    FVTMd = wgms3d_int(abs(FVTM.hr).^2+abs(FVTM.hz).^2+abs(FVTM.hp).^2, FVTM.r, FVTM.z);
    dTEhr(lambdac) = wgms3d_int((abs(FVTE.hr)-abs(SVTE.hr)).^2, FVTE.r, FVTE.z) / FVTEd;
    dTEhz(lambdac) = wgms3d_int((abs(FVTE.hz)-abs(SVTE.hz)).^2, FVTE.r, FVTE.z) / FVTEd;
    dTEhp(lambdac) = wgms3d_int((abs(FVTE.hp)-abs(SVTE.hp)).^2, FVTE.r, FVTE.z) / FVTEd;
    dTMhr(lambdac) = wgms3d_int((abs(FVTM.hr)-abs(SVTM.hr)).^2, FVTM.r, FVTM.z) / FVTMd;
    dTMhz(lambdac) = wgms3d_int((abs(FVTM.hz)-abs(SVTM.hz)).^2, FVTM.r, FVTM.z) / FVTMd;
    dTMhp(lambdac) = wgms3d_int((abs(FVTM.hp)-abs(SVTM.hp)).^2, FVTM.r, FVTM.z) / FVTMd;

    semivectorial_plot
    drawnow
    
end
    

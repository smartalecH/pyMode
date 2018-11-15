
% It is well-known that some waveguide structures, even when they consist
% only of lossless dielectrics and perfectly conducting metals, may have
% modes with a complex (i.e., neither purely real nor purely imaginary)
% propagation constant. This can happen even if the waveguide is straight
% and does not have any obvious means for energy to leak away from the
% waveguide core.

% wgms3d can also be used to calculate these 'complex modes'. A PML is not
% necessary. Just run the mode solver with the respective geometry; if
% complex modes are there, they will be returned just like any other mode.

% Here, as an example, we consider the waveguide from Fig. 4b in J. Strube,
% F. Arndt, "Rigorous Hybrid-Mode Analysis of the Transition from
% Rectangular Waveguide to Shielded Dielectric Image Guide," IEEE
% Transactions on Microwave Theory and Techniques, vol. MTT-33, no. 5, May
% 1985. This script reproduces the data from Fig. 4b.

fs = linspace(12e9, 18e9, 200);

ms = 'LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -d -U -7.899e-3:250:7.899e-3 -V 0:250:7.9e-3 -n 5';

neffs = nan(5, length(fs));

figure(1); clf
drawnow

for fc = 1 : length(fs)
        
    lambda = 3e8 / fs(fc);

    com = sprintf('%s -l %e -g mgp.mgp', ms, lambda);
    m = wgms3d_run(com);

    for i = 1 : size(neffs,1)
        neffs(i,fc) = m{i}.neff;
    end

    strube_plot
    drawnow
    
end
    

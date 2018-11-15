
% Reproduces results from Fig. 4 of this article:
%
% Thach G. Nguyen, Ravi S. Tummidi, Thomas L. Koch, Arnan Mitchell,
% "Rigorous Modeling of Lateral Leakage Loss in SOI Thin-Ridge Waveguides
% and Couplers", IEEE Photonics Technology Letters, vol. 21, no. 7, April
% 2009.

nsearch = [ 1.721 ];

tic

[ Ws, Modes ] = wgms3d_tracemodes (1.36e-6, nsearch, ...
                                   0.002e-6, 0.0001e-6, 1.48e-6, ...
                                   @nguyen_prepfunc, @nguyen_plotfunc);

disp(sprintf('Calculation took %ds = %dmin = %dh.', ...
	     ceil(toc), ceil(toc/60), ceil(toc/60/60)))

wgms3d_clean


lambda0 = 1.55;

% Here are some different discretization grids:

% bendit01.mat:
nsearch = [ 1.4645060481510088 1.4526985903868901 ];
ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -U -6:200:4.1 -V -5:200:2.51 -P e:30:0.4 -P s:30:1.0';

% bendit02.mat:
% nsearch = [ 1.4645 1.4526 ];
% ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -U -6:200:8 -V -5:300:+2.5 -P e:30:1.0 -P s:30:1.0';

% bendit03.mat:
% nsearch = [ 1.4645 1.4526 ];
% ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -U -4:400:6 -V -3.9:300:+2.5 -P e:30:0.7 -P s:30:1.0';

% bendit04.mat:
% nsearch = [ 1.4645 1.4526 ];
% ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -U -5.1:650:6 -V -3.951:400:+2.5 -P e:27:0.8 -P s:30:0.9';

% bendit05.mat:
% nsearch = [ 1.4645 1.4526 ];
% ms = 'LD_LIBRARY_PATH=; ~/wgms3d/wgms3d-0.6/wgms3d -U -5.05:601:6 -V -3.971:413:+2.477 -P e:28:0.67 -P s:30:0.95';

% bendit06.mat:
% nsearch = [ 1.4645 1.4526 ];
% ms = 'LD_LIBRARY_PATH=; ~/wgms3d/wgms3d-0.6/wgms3d -U -5.0:500:5.7 -V -3.9:401:+2.4 -P e:15:0.8 -P s:20:0.8';

% bendit07.mat:
%nsearch = [ 1.4645 1.4526 ];
%ms = 'LD_LIBRARY_PATH=; ~/wgms3d/wgms3d-0.6/wgms3d -U -5.0:500:6.0 -V -4.0:400:+2.5 -P e:15:0.8 -P s:15:0.8';

%bendit08.mat:
%nsearch = [ 1.4645 1.4526 ];
%ms = 'LD_LIBRARY_PATH=; ~/wgms3d/wgms3d-0.6/wgms3d -U -4.5:600:6.5 -V -4.5:450:+2.5 -P e:20:1.0 -P s:20:1.0';

% bendit09.mat:
%nsearch = [ 1.4645 1.4526 ];
%ms = 'LD_LIBRARY_PATH=; ~/wgms3d/wgms3d-0.6/wgms3d -U -5.0:550:6.5 -V -4.0:500:+3.0 -P e:20:0.5 -P s:20:0.5';

% bendit09b.mat:
%nsearch = [ 1.4645 1.4526 ];
%ms = 'LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -U -5.0:550:6.5 -V -4.0:500:+3.0 -P e:20:1.0 -P s:20:1.0';

prepfunc = @(Rcurv,nn) (sprintf('%s -g prkna.mgp -n 4 -E -F -l %e -s %e -R %e', ...
                                ms, lambda0, abs(mean(nn)), Rcurv));
plotfunc = @bendit_plotfunc;

tic

[ Rcurvs, Modes ] = wgms3d_tracemodes (50.0, nsearch, ...
                                       -5.0, 0.1, 10.0, ...
                                       prepfunc, plotfunc);

disp(sprintf('Calculation took %ds = %dmin = %dh.', ...
	     ceil(toc), ceil(toc/60), ceil(toc/60/60)))

wgms3d_clean

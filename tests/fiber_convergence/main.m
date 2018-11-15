
ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -l 10 -g fiber.mgp';
symm = '-n 2'
neffreal = [ 2.161722239740193e+00
             2.102237229347229e+00 ];
calc
save case13_1_dnonzero.mat

ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -l 10 -g fiber.mgp';
symm = '-n 2 -M s'
neffreal = [ 2.351757217372069e+00
             1.843555913442746e+00 ];
calc
save case13_2_dnonzero.mat

ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -l 10 -g fiber.mgp';
symm = '-n 2 -M w'
neffreal = [ 2.351757217372069e+00
             1.843555913442746e+00 ];
calc
save case13_3_dnonzero.mat

ms = 'LD_LIBRARY_PATH=; ~/wgms3d/current/wgms3d -l 10 -g fiber.mgp';
symm = '-n 4 -M w -M s'
neffreal = [ 2.102237229347229e+00
             2.082244850080460e+00 ];
calc
save case13_4_dnonzero.mat

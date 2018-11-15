
% Keep CD constant, vary grid size

Ns = 150;
ms = 'LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -g t02_fiberum.mgp';

neffreal = 1.979061369181741e+00;

modes = wgms3d_run(sprintf(['%s -l 1.55 -e -U 0:%d:1.5 -V 0:%d:1.5 -n 1 -E ' ...
                    '-F -G -H -M s'], ms, Ns, Ns));
modes{1} = wgms3d_load_mode_field(modes{1});

x=real(modes{1}.r);
y=real(modes{1}.z);

if 1
    cfig(1); clf
    subplot(2,1,1)
    surf(x,y,sqrt(modes{1}.hr.^2 + modes{1}.hz.^2));
    shading interp;
    title 'numerical'
    subplot(2,1,2)
    surf(x,y,imag(modes{1}.hp))
    shading interp;
end

if 0
    [exr, eyr, ezr] = circular_exact_e_field (1, neffreal, x, y, .160, 1.55, 3.5, 1.5);
    fact = modes{1}.er(1,1) / exr(1,1);
end
if 1
    [exr, eyr, ezr] = circular_exact_h_field (1, neffreal, x, y, .160, 1.55, 3.5, 1.5);
    fact = modes{1}.hz(1,1) / eyr(1,1);
end
exr = exr * fact;
eyr = eyr * fact;
ezr = imag(ezr * fact);
if 1
    cfig(2); clf
    subplot(2,1,1)
    surf(x,y,sqrt(exr.^2 + eyr.^2));
    shading interp;
    title 'exact'
    subplot(2,1,2)
    surf(x,y,ezr)
    shading interp;
end

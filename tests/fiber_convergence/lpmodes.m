
% (früher: w10.m)

% Keep CD constant, vary grid size

Ns = linspace(50,400,25);
Ns = floor(Ns);

ms = 'LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -d -l 10 -p -g fiber.mgp';

neffreal = [ 2.370506655212511e+00
             2.161722239629984e+00
             2.161722239629984e+00
             1.865029549757501e+00
             1.865029549757501e+00
             1.770543445693141e+00 ];

ls = { 'x-', 'xr-', 'xr-', 'x--', 'x--', 'xk-' };
  
neffs = nan * ones(length(neffreal), length(Ns));

for Nc = 1 : length(Ns)
  com = sprintf('%s -U -9.99:%d:9.99 -V -9.99:%d:9.99 -n 6', ms, Ns(Nc), Ns(Nc));
  modes = wgms3d_run(com);
  
  for i = 1 : 6
      neffs(i,Nc) = modes{i}.neff;
  end
  
  cfig(4); clf

  for i = 1 : 6
    relerrs = abs(neffs(i,:) - neffreal(i)) / neffreal(i);
    loglog(1 ./ Ns, relerrs, ls{i})
    hold on
  end
  xlabel 'Resolution (inverse grid-point distance) [a.u.]'
  ylabel 'Relative error in n_{eff}'
  title(strrep(sprintf('%s on %s', mfilename('fullpath'), date), '_', '\_'))
  ylim([10^-8 1])
  
  drawnow
end

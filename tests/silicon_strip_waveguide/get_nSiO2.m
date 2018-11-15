
% Wavelength = 'l' in microns

function n = get_nSiO2 (l)
  
  BSiO21 = 0.6961663;
  BSiO22 = 0.4079426;
  BSiO23 = 0.8974794;
  CSiO21 = 0.004679148;
  CSiO22 = 0.01351206;
  CSiO23 = 97.934002;

  % Das Buch, S. 6 -> von 0.21µ bis 3.71µ .. ok!
  n = sqrt(1 + BSiO21*l^2/(l^2-CSiO21) + BSiO22*l^2/(l^2-CSiO22) + ...
           BSiO23*l^2/(l^2-CSiO23));

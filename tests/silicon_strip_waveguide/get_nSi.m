
% Wavelength = 'l' in microns

function n = get_nSi (l)
  
  eSi1 = 11.6858;
  ASi = 0.939816;
  BSi = 8.10461e-3;
  lSi1 = 1.1071;

  % Dimi OPEX 2004
  % Sellmeier-Formel von "Li" -> 1.2 bis 14µ .. ok!
  n = sqrt(eSi1 + ASi/(l^2) + BSi*l^2/(l^2-lSi1^2));


function [ex, ey, ez] = circular_exact_e_field (nu, neff, xs, ys, ...
                                                rho, lambda, nco, ncl)
  
  k = 2*pi/lambda;
  V = k*rho*sqrt(nco^2-ncl^2);
  beta = k * neff;
  U = rho * sqrt(k^2*nco^2 - beta.^2);
  W = rho * sqrt(beta.^2 - k^2*ncl^2);
  I = sqrt(-1);
  
  u0 = 1.256e-6;
  e0 = 8.854e-12;

  b1 = 1/(2*U)*(besselj(nu-1,U)/besselj(nu,U) - besselj(nu+1,U)/besselj(nu,U));
  b2 = -1/(2*W)*(besselk(nu-1,W)/besselk(nu,W) + besselk(nu+1,W)/besselk(nu,W));
  delta = (nco^2-ncl^2)/(2*nco^2);
  F1 = (U*W/V)^2*(b1+(1-2*delta)*b2)/nu;
  F2 = (V/(U*W))^2*nu/(b1+b2);
  a1 = (F2-1)/2;
  a2 = (F2+1)/2;
  a3 = (F1-1)/2;
  a4 = (F1+1)/2;
  a5 = (F1-1+2*delta)/2;
  a6 = (F1+1-2*delta)/2;
  
  for j = 1 : length(ys)
    y = ys(j);
    for i = 1 : length(xs)
      x = xs(i);
      r = sqrt(x^2 + y^2);
      phi = atan2(y, x);
      R = r/rho;
      
      if R < 1
        er = -(a1*besselj(nu-1,U*R) + a2*besselj(nu+1,U*R))/besselj(nu,U) ...
             * cos(nu*phi);
        ephi = -(a1*besselj(nu-1,U*R) - a2*besselj(nu+1,U*R))/besselj(nu,U) ...
               * (-sin(nu*phi));
        ez(j,i) = -I*U/(rho*beta) * besselj(nu,U*R) / besselj(nu,U) * cos(nu*phi);
      else
        er = -(U/W)*(a1*besselk(nu-1,W*R) - a2*besselk(nu+1,W*R))/ ...
             besselk(nu,W) * cos(nu*phi);
        ephi = -(U/W)*(a1*besselk(nu-1,W*R) + a2*besselk(nu+1,W*R))/ ...
               besselk(nu,W) * (-sin(nu*phi));
        ez(j,i) = -I*U/(rho*beta) * besselk(nu,W*R) / besselk(nu,W) * cos(nu*phi);
      end
      
      ex(j,i) = er*cos(phi) + ephi*(-sin(phi));
      ey(j,i) = er*sin(phi) + ephi*cos(phi);
    end
  end


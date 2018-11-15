
%   wgms3d - a full-vectorial finite-difference mode solver.
%
%   Copyright (C) 2005-2010  Michael Krause <m.krause@tu-harburg.de>
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

function wgms3d_plot_mgp (fn)

    f = fopen(fn, 'r');
    if f <= 0
        f = fopen([fn '.mgp'], 'r');
        if f <= 0
            error('Can''t read file')
        end
    end

    while ~feof(f)
        
        s = fgets(f);
        
        switch s(1)
          case 'x'
            break

          case 'l'
            a = sscanf(s, 'l %*s %*s %lf %lf %lf %lf');
            plot([a(1) a(3)], [a(2) a(4)], 'Color', [ 0.8 0.3 0.0 ])
            hold on
            continue
            
          case 'b'
            a = sscanf(s, 'b %*s %*s %lf %lf %lf %lf %lf %lf');
            u = linspace(0, 1, 20);
            if 0
                Bx = (1-u).^2*a(1) + 2*u.*(1-u)*a(3) + u.^2*a(5);
                By = (1-u).^2*a(2) + 2*u.*(1-u)*a(4) + u.^2*a(6);
                plot(Bx, By, 'Color', [ 0.8 0.3 0.0 ])
                hold on
            end
            if 1
                Bx = (1-u).^3.*(1+u)*a(1) + 2*u.*(1-u).*(-u.^2+u+1)*a(3) + u.^3.*(2-u)*a(5);
                By = (1-u).^3.*(1+u)*a(2) + 2*u.*(1-u).*(-u.^2+u+1)*a(4) + u.^3.*(2-u)*a(6);
                plot(Bx, By, 'Color', [ 0.8 0.3 0.0 ])
                hold on
            end
            continue
            
          case 'c'
            a = sscanf(s, 'c %*s %*s %lf %lf %lf');
            u = linspace(0, 2*pi, 51);
            Bx = a(1) + a(3)*cos(u);
            By = a(2) + a(3)*sin(u);
            plot(Bx, By, 'Color', [ 0.8 0.3 0.0 ])
            hold on
            continue
            
          case 'e'
            a = sscanf(s, 'e %*s %*s %lf %lf %lf %lf');
            u = linspace(0, 2*pi, 51);
            Bx = a(1) + a(3)*cos(u);
            By = a(2) + a(4)*sin(u);
            plot(Bx, By, 'Color', [ 0.8 0.3 0.0 ])
            hold on
            continue
            
          case 'n'
          case 'y'
          case 'C'
          case '#'
          case '%'
            continue
            
          otherwise
            error 'Can''t parse MGP file'
            
        end
    end
    
    fclose(f);
    
    hold off

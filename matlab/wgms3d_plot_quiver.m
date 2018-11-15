
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

% Plots the specified vector function ux,uy(x,y) as a quiver plot.

function wgms3d_plot_quiver (x0, y0, ex, ey, varargin)
    
    maxabsreal = max(max(sqrt(real(ex).^2+real(ey).^2)));
    maxabsimag = max(max(sqrt(imag(ex).^2+imag(ey).^2)));
    if maxabsimag <= 1e-3 * maxabsreal
        % Plot only real part:
        wgms3d_plot_quiver_sub(x0, y0, real(ex), real(ey), 1, [0 0 1], varargin{:})
    elseif maxabsreal <= 1e-3 * maxabsimag
        % Plot only imaginary part:
        wgms3d_plot_quiver_sub(x0, y0, real(ex), real(ey), 1e-6, [0 0 1], varargin{:})
        wgms3d_plot_quiver_sub(x0, y0, imag(ex), imag(ey), 1, [0 0.7 0], varargin{:})
    else
        % Plot both, correctly scaled:
        if maxabsreal > maxabsimag
            wgms3d_plot_quiver_sub(x0, y0, real(ex), real(ey), 1, [0 0 1], varargin{:})
            wgms3d_plot_quiver_sub(x0, y0, imag(ex), imag(ey), maxabsimag/maxabsreal, [0 0.7 0], varargin{:})
        else
            wgms3d_plot_quiver_sub(x0, y0, real(ex), real(ey), maxabsreal/maxabsimag, [0 0 1], varargin{:})
            wgms3d_plot_quiver_sub(x0, y0, imag(ex), imag(ey), 1, [0 0.7 0], varargin{:})
        end
    end
end

function wgms3d_plot_quiver_sub (x, y, ux, uy, scale, color, varargin)

    wgms3d_plot_parse_args

    if ~doquiver
        return
    end
  
    if quiver_normed
        ur = sqrt(real(ux).^2+real(uy).^2);
        ux = real(ux) ./ ur;
        uy = real(uy) ./ ur;
        if scale > 1e-5
            scale = 1;
        end
    end
  
    if xsymm > 0
        x = [ x ; -flipud(x) ];
        if xor(xsymm == 2, 1)
            ux = [ ux fliplr(ux) ];
            uy = [ uy -fliplr(uy) ];
        else
            ux = [ ux -fliplr(ux) ];
            uy = [ uy fliplr(uy) ];
        end
    end

    if quiver_grid > 0
        nx = quiver_grid;
        ny = round((y(end)-y(1)) / (x(end)-x(1)) * nx);
        xo = x;
        yo = y;
        x = linspace(xo(1), xo(end), nx)';
        y = linspace(yo(1), yo(end), ny)';
        [Xo,Yo] = meshgrid(xo,yo);
        ux = interp2(Xo, Yo, ux, x, y');
        uy = interp2(Xo, Yo, uy, x, y');
    end
  
    hold on
    quiver(x, y, real(ux), real(uy), scale * quiver_scale, 'LineWidth', quiver_lw, 'Color', color)
end

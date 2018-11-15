
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

function wgms3d_plot_post (xcomplex, ycomplex, varargin)
  
    wgms3d_plot_parse_args
  
    x0 = real(xcomplex);
    y0 = real(ycomplex);
    
    if xsymm
        x0 = [ x0 ; -flipud(x0) ];
    end

    % Plot geometry
    if length(mgp) > 0
        hold on
        wgms3d_plot_mgp(mgp)
    end

    % Plot PML boundaries
    pmls_x = wgms3d_find_pml(imag(xcomplex));
    pmls_y = wgms3d_find_pml(imag(ycomplex));

    hold on
    for i = 1 : length(pmls_x)
        plot([x0(pmls_x(i)) x0(pmls_x(i))], [y0(1) y0(end)], ...
             'k--', 'LineWidth', 3)
    end
    for i = 1 : length(pmls_y)
        plot([x0(1) x0(end)], [y0(pmls_y(i)) y0(pmls_y(i))], ...
             'k--', 'LineWidth', 3)
    end

    % Plot CD boundary
    wgms3d_plot_cd_boundary(x0, y0, [0.8 0.0 0.0])

    hold off
    

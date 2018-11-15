
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

function wgms3d_plot_ht (mode, varargin)

    wgms3d_plot_parse_args
  
    data = { 'hr', 'hz' };
    wgms3d_plot_get_data

    if realcontours
        forcontours = sqrt(real(hr).^2+real(hz).^2);   
    else
        forcontours = sqrt(abs(hr).^2 + abs(hz).^2);
    end

    clf
    wgms3d_plot_contours(x0, y0, forcontours, varargin{:})
    wgms3d_plot_quiver(x0, y0, hr, hz, varargin{:})
    wgms3d_plot_post(xcomplex, ycomplex, varargin{:})
 
    title(sprintf('Transverse magnetic field of mode %d', cnt))
    xlabel '\rho'
    ylabel 'z'

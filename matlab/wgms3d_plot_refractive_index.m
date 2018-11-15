
%   wgms3d - a full-vectorial finite-difference mode solver.
%
%   Copyright (C) 2005-2011  Michael Krause <m.krause@tu-harburg.de>
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

function wgms3d_plot_refractive_index

    x = load('r.txt');
    y = load('z.txt');
    f = fopen('epsis.bin');
    epsis = wgms3d_read_bin_data(x, y, f);
    fclose(f);

    x = x(:,1);
    y = y(:,1);
    n = sqrt(epsis);

    clf
    surface(x, y, abs(n), 'EdgeColor', 'none')
    colormap gray

    title('Refractive-index profile')
    xlabel '\rho'
    ylabel 'z'

    hold on
    wgms3d_plot_cd_boundary(x, y, [0.8 0 0])
    hold off
    
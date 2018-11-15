
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

function wgms3d_plot_cd_boundary (x0, y0, color)

    plot([x0(1) x0(end) x0(end) x0(1) x0(1)], ...
         [y0(1) y0(1) y0(end) y0(end) y0(1)], ...
         '-', 'LineWidth', 3, 'Color', color)
    
    x1 = x0(1);
    x2 = x0(end);
    y1 = y0(1);
    y2 = y0(end);

    wx = x2-x1;
    wy = y2-y1;
    cx = (x1+x2)/2;
    cy = (y1+y2)/2;
  
    xlim([cx-0.55*wx cx+0.55*wx])
    ylim([cy-0.55*wy cy+0.55*wy])

    
    
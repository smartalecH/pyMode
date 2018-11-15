
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


% This function calculates the group index of a mode from its field
% distribution at a single wavelength. Material dispersion is ignored. See
% Snyder/Love, Eq. (31-31).

% Example usage:
% >> m = wgms3d_run('LD_LIBRARY_PATH= ~/wgms3d/current/wgms3d -g fiber.mgp -l .1 -n 5 -U -1:200:+1 -V -1:200:+1 -E -F -G -e');
% >> m{1} = wgms3d_load_mode_field(m{1});
% >> wgms3d_group_index_from_field(m{1})

function ng = wgms3d_calc_group_index (mode)

    eps0 = 8.85418781762e-12;
    c0 = 2.99792458e8;
    
    N = wgms3d_modeproduct(mode, mode);
    absEsq = abs(mode.er).^2 + abs(mode.ez).^2 + abs(mode.ep).^2;
    D = (eps0 / 2) * wgms3d_int(mode.epsis .* absEsq, mode.r, mode.z);
    vg = N / D;
    
    ng = c0 / vg;

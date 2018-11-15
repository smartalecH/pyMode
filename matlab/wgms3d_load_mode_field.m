
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

function mode = wgms3d_load_mode_field (mode)

    dir = '.';
    
    rfn = sprintf('%s/r.txt', dir);
    r = load(rfn);
    if isempty(r)
        error(sprintf('''%s'' not found', rfn))
    end
    mode.r = r(:,1) + i*r(:,2);
  
    zfn = sprintf('%s/z.txt', dir);
    z = load(zfn);
    if isempty(z)
        error(sprintf('''%s'' not found', zfn))
    end
    mode.z = z(:,1) + i*z(:,2);

    try
        f = fopen(sprintf('%s/epsis.bin', dir), 'r');
        mode.epsis = wgms3d_read_bin_data(mode.r, mode.z, f);
        fclose(f);
    end
    
    try
        [ mode.er, mode.beta ] = wgms3d_load_field (mode.r, mode.z, ...
                                        sprintf('%s/er', dir), mode.nev);
    end

    try
        [ mode.ez, mode.beta ] = wgms3d_load_field (mode.r, mode.z, ...
                                        sprintf('%s/ez', dir), mode.nev);
    end

    try
        [ mode.ep, mode.beta ] = wgms3d_load_field (mode.r, mode.z, ...
                                        sprintf('%s/ep', dir), mode.nev);
    end

    try
        [ mode.hr, mode.beta ] = wgms3d_load_field (mode.r, mode.z, ...
                                        sprintf('%s/hr', dir), mode.nev);
    end

    try
        [ mode.hz, mode.beta ] = wgms3d_load_field (mode.r, mode.z, ...
                                        sprintf('%s/hz', dir), mode.nev);
    end

    try
        [ mode.hp, mode.beta ] = wgms3d_load_field (mode.r, mode.z, ...
                                        sprintf('%s/hp', dir), mode.nev);
    end

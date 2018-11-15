
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

function [ data ] = wgms3d_read_bin_data (x, y, f)
    
    data = fread(f, [ 2*length(x) length(y) ], 'double');
    real_part_indices = 1:2:2*length(x);
    imag_part_indices = real_part_indices + 1;
    data = data(real_part_indices,:)' + i*data(imag_part_indices,:)';

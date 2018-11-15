
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

% Make the mode-field components listed in cell-array 'data' available in
% the current scope. Can handle modes passed as a structure to the plot
% function or modes not loaded into workspace but residing on the
% filesystem.

if class(mode) == 'struct'
    xcomplex = mode.r;
    ycomplex = mode.z;
    x0 = real(xcomplex);
    y0 = real(ycomplex);
    
    for j = 1 : length(data)
        eval(sprintf('%s = mode.%s;', data{j}, data{j}));
    end
    
    cnt = mode.nev;
else
    xcomplex = load('r.txt');
    ycomplex = load('z.txt');
    xcomplex = xcomplex(:,1) + sqrt(-1)*xcomplex(:,2);
    ycomplex = ycomplex(:,1) + sqrt(-1)*ycomplex(:,2);
    x0 = real(xcomplex);
    y0 = real(ycomplex);
        
    for j = 1 : length(data)
        eval(sprintf('[ %s, beta ] = wgms3d_load_field(x0, y0, ''%s'', mode);', ...
                     data{j}, data{j}));
    end

    cnt = mode;
end

if pphase ~= 0
    for j = 1 : length(data)
        eval(sprintf('%s = %s .* exp(sqrt(-1) * pphase);', data{j}, data{j}));
    end
end

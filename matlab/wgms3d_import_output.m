
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


%   Function to convert wgms3d's shell output into Matlab variables.


function modes = wgms3d_import_output (r)
    
    modes0 = {};
    
    n = 1;
    abbrechen = 0;

    while true

        clear mode

        while true
            [l,r] = wgms3d_sgetl(r);
            if isempty(l)
                abbrechen = 1;
                break
            end
            if ~isempty(strfind(l, 'EV'))
                a = sscanf(l, 'EV %d: n_eff = %f + i%e');
                mode.nev = a(1);
                mode.neff = a(2) + sqrt(-1)*a(3);
                break
            end
        end

        if abbrechen == 1
            break
        end
        
        neffs(n) = mode.neff;

        [l,r] = wgms3d_sgetl(r);
        b = sscanf(l, '        alpha = %fdB/UOL [%fdB/90deg], pol = ''%c''');
        mode.alphadbm = b(1);
        mode.alphadb90 = b(2);
        mode.pol = b(3);
        
        modes0{n} = mode;
        n = n + 1;
    end
    
    % Now sort modes0 by descending neff^2
    [a,b] = sort(real(neffs.^2), 2, 'descend');
    for i = 1 : length(modes0)
        modes{i} = modes0{b(i)};
    end

    
function [l,s] = wgms3d_sgetl (s)

    p = findstr(s, 10);
  
    if ~isempty(p)
        l = s(1:p(1));
        s = s(p(1)+1:end);
    else
        l = [];
        s = [];
    end

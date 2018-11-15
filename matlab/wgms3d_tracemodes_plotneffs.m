
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

function wgms3d_tracemodes_plotneffs (PARs, Modes)
    
    for i = 1 : length(PARs)
        for j = 1 : length(Modes{i})
            neffs(j,i) = abs(real(Modes{i}{j}.neff));
            alphas(j,i) = abs(Modes{i}{j}.alphadbm);
        end
    end
    
    ccfig(9); clf
    subplot(2,1,1)
    plot(PARs, neffs, 'x-')
    ylabel 'Re n_{eff}'
    xlabel 'PAR'
    subplot(2,1,2)
    semilogy(PARs, alphas, 'x-')
    ylabel '\alpha [dB/UOL]'
    xlabel 'PAR'

    return

function ccfig (h)
  
    if ishandle(h)
        set(0,'CurrentFigure',h)
    else
        figure(h)
    end


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

if PAR == PARdest
    % We're done!
    finished = true;
    return;
end

finished = false;
        
% Update tracing parameter
PAR = PARs(end) + PARstep;
if (PARstep > 0 && PAR > PARdest) ...
        || (PARstep < 0 && PAR < PARdest) ...
        || abs(PAR - PARdest) < abs(PARstep_min)/100
    PAR = PARdest;
    PARstep = PARdest - PARs(end);
end

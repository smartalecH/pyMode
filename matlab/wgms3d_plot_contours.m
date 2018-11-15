
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

% Plots the specified scalar positive function u(x,y) either as
% contours or as a 3-D plot, as desired.

function wgms3d_plot_contours (x, y, u, varargin)
  
    wgms3d_plot_parse_args
  
    if xsymm > 0
        x = [ x ; -flipud(x) ];
        u = [ u fliplr(u) ];
    end
  
    if ~threedee
        %mmin = min(min(u));
        mmax = max(max(u));
        if logcontours
            % 1dB-Contours (1dB = Faktor .7943)
            % 2dB-Contours (2dB = Faktor .631)
            %      max(max(u))
            clevels = logspace(log10(mmax), ...
                               log10(mmax * (10^(-logcontours/10))^(ncontours-1)), ...
                               ncontours + 1);
            clevels = clevels(2:end);
        else
            clevels = linspace(mmax, mmax/ncontours, ncontours + 1);
        end
    
        if 0
            % old version:
            [C,h,CF] = contourf(x, y, u, clevels);
        else
            % for Matlab 7 and later:
            [C,h] = contourf(x, y, u, clevels);
            h = get(h,'Children');
        end
        %    colormap 'hot';
        colormap 'gray';
        X = (length(colormap)-1)/(length(clevels)-1);
        for i = 1 : length(h)
            ud = get(h(i),'UserData');
	    if ~isempty(ud)
		% this will only work in Matlab:
		a = find(clevels == ud);
		ci = floor(1 + (a-1)*X);
		set(h(i), 'CData', ci)
		set(h(i), 'FaceVertexCData', ci)
		set(h(i), 'CDataMapping', 'direct')
		set(h(i), 'EdgeColor', edgecolor)
	    end
        end
        
    else
        surf(x, y, u);
        shading interp;
        colormap([.8 .8 .8])
        lighting phong
        light('Position', [-7 7 3])
    end
    
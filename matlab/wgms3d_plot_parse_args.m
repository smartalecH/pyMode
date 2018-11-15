
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

xsymm = false;
mgp = '';
doquiver = true;
swapxy = false;
quiver_grid = 80;
quiver_normed = false;
quiver_scale  = 1.0;
quiver_lw = 1.0;
ncontours = 10;
threedee = false;
logcontours = 0;
realcontours = false;
edgecolor = 'none';
pphase = 0;
mmax = 0;

i = 1;
while i <= length(varargin)
  if strcmp(varargin(i), 'Symmetry')
    if strcmp(varargin(i+1), 'X')
      xsymm = 1;
    elseif strcmp(varargin(i+1), 'x')
      xsymm = 2;
    else
      xsymm = false;
    end
    i = i + 2;
    continue
  end
  if strcmp(varargin(i), 'Mode')
    if strcmp(varargin(i+1), 'Flat')
      threedee = false;
    elseif strcmp(varargin(i+1), '3D')
      threedee = true;
    else
      error 'Can''t handle argument list.'
    end
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'Geometry')
    mgp = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'NoQuiver')
    i = i + 1;
    doquiver = false;
    continue;
  end
  if strcmp(varargin(i), 'QuiverNormed')
    i = i + 1;
    quiver_normed = true;
    continue;
  end
  if strcmp(varargin(i), 'QuiverScale')
    quiver_scale = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'QuiverGrid')
    quiver_grid = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'QuiverLineWidth')
    quiver_lw = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'Contours')
    ncontours = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'RealContours')
    i = i + 1;
    realcontours = true;
    continue;
  end
  if strcmp(varargin(i), 'LogContours')
    logcontours = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'EdgeColor')
    edgecolor = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'Phase')
    pphase = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  if strcmp(varargin(i), 'Max')
    mmax = cell2mat(varargin(i+1));
    i = i + 2;
    continue;
  end
  
  error(sprintf('Can''t handle argument list: %s', cell2mat(varargin(i))))
end

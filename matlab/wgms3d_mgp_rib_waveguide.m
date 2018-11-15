
%   wgms3d - a full-vectorial finite-difference mode solver.
%
%   Copyright (C) 2005-2012  Michael Krause <m.krause@tu-harburg.de>
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

% Generates a geometry like this:
%
%              /-----------\
%              |           |
%   -----------/           \---------------------
%   ---------------------------------------------
%
% w = central rib width
% h = central rib height
% d = etch depth (thus, the height of the adjacent slabs is given by 'h-d')
%
% Setting d=h in order to get a completely etched-down waveguide works, but
% corners will be rounded like this:
%
%    /--------\
%    |        |
% ___/________\____
%
% instead of
%
%    /--------\
%    |        |
% ___\________/____
%
% which might be more what you want. In that case, specify an arbitrary
% argument in addition to nBuf.
%
% For example usage, see tests/2006_OPEX_Dulkeith/ and tests/lossy_materials

function wgms3d_mgp_rib_waveguide (fn, w, h, d, rc, nTop, nCore, nBuf, varargin)
  
    assert(w > 0);
    assert(d > 0);
    assert(h > 0);
    assert(d <= h, 'Etch depth larger than rib height.');
    assert(d >= 2*rc, 'Etch depth not large enough to leave room for curved corners.');
    
    assert(nCore > nBuf);
    assert(nCore > nTop);
    
    s = h - d; % height of the two slabs adjacent to the central rib
    xmax = 100*w;

    % Round two lower corners outside or inside?
    if length(varargin) == 0
        xc = w/2 + rc;
    else
        assert(d == h);
        xc = w/2 - rc;
    end

    f = fopen(fn, 'w');

    fprintf(f, 'n (%e,%e)\n', real(nTop), imag(nTop));

    % Draw upper three straight boundaries (core <-> top-cladding)
    fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
            real(nTop), imag(nTop), real(nCore), imag(nCore), -w/2, s+rc, -w/2, h-rc);
    fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
            real(nTop), imag(nTop), real(nCore), imag(nCore), -w/2+rc, h, +w/2-rc, h);
    fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
            real(nTop), imag(nTop), real(nCore), imag(nCore), +w/2, h-rc, +w/2, s+rc);

    % Now draw lower straight boundary and top-cladding-buffer-cladding interface
    if s > 0
        % Rib waveguide
        fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
                real(nTop), imag(nTop), real(nCore), imag(nCore), -xmax, s, -xc, s);
        fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
                real(nTop), imag(nTop), real(nCore), imag(nCore), +xc, s, +xmax, s);
        fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
                real(nCore), imag(nCore), real(nBuf), imag(nBuf), -xmax, 0, +xmax, 0);
    else
        % Strip waveguide
        if nTop ~= nBuf
            fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
                    real(nTop), imag(nTop), real(nBuf), imag(nBuf), -xmax, 0, -xc, 0);
            fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
                    real(nTop), imag(nTop), real(nBuf), imag(nBuf), +xc, 0, +xmax, 0);
        end
        fprintf(f, 'l (%e,%e) (%e,%e) %e %e %e %e\n', ...
                real(nCore), imag(nCore), real(nBuf), imag(nBuf), -xc, 0, +xc, 0);
    end

    % Now draw rounded corners, if requested
    if rc > 0
        % Upper two rounded corners
        fprintf(f, 'b (%e,%e) (%e,%e) %e %e %e %e %e %e\n', ...
                real(nTop), imag(nTop), real(nCore), imag(nCore), -w/2, h-rc, -w/2, h, -w/2+rc, h);
        fprintf(f, 'b (%e,%e) (%e,%e) %e %e %e %e %e %e\n', ...
                real(nTop), imag(nTop), real(nCore), imag(nCore), +w/2-rc, h, +w/2, h, +w/2, h-rc);
        % Lower two rounded corners
        fprintf(f, 'b (%e,%e) (%e,%e) %e %e %e %e %e %e\n', ...
                real(nTop), imag(nTop), real(nCore), imag(nCore), -xc, s, -w/2, s, -w/2, s+rc);
        fprintf(f, 'b (%e,%e) (%e,%e) %e %e %e %e %e %e\n', ...
                real(nTop), imag(nTop), real(nCore), imag(nCore), +w/2, s+rc, +w/2, s, +xc, s);
    end
    
    fprintf(f, 'C %e %e %e %e\n', ...
            -w/2, 0, +w/2, h);

    fprintf(f, 'x\n');
    fclose(f);

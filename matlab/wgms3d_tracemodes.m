
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

function [ PARs, Modes ] = wgms3d_tracemodes (PAR, ini, ...
                                              PARstep_orig, PARstep_min, PARdest, ...
                                              prepfunc, plotfunc)

    PARstep = PARstep_orig;
    
    if iscell(ini)
        % User has specified a set of modes including mode-field data which are to
        % be used as a starting point. Those modes are again returned as the
        % first entry in the results.
        Modes = ini;
        PARs = PAR;
        for i = 1 : length(Modes{1})
            nsearch(i) = abs(real(Modes{1}{i}.neff));
        end
        wgms3d_tracemodes_update_parameter;
        if finished
            return
        end
    else
        % User has specified a set of effective indices to be used as a starting
        % point.
        Modes = {};
        PARs = [];
        nsearch = ini;
    end
    
    stepok = 0;
    
    while 1

        steplowered = false;

        wgms3d_clean
        % Call prepfunc. This function can update the geometry file based on the new
        % parameter 'PAR', it can generate new discretization meshes, and it
        % returns the new wgms3d command line.
        com = prepfunc(PAR, nsearch);
        m = wgms3d_run(com);

        if ~isempty(m)
            
            for i = 1 : length(m)
                m{i} = wgms3d_load_mode_field(m{i});
                m{i}.K = wgms3d_modeproduct(m{i},m{i});
            end

            clear NewModes;
            verbastelt = [];
            modelost = false;
            
            if isempty(Modes)

                % First run, no mode fields to compare against
                for mc = 1 : length(nsearch)
                    ds = [];
                    for i = 1 : length(m)
                        ds(i) = abs(abs(m{i}.neff) - nsearch(mc));
                        if length(find(verbastelt == i)) == 1
                            ds(i) = inf;
                        end
                    end
                    ds
                    a = find(ds == min(ds)); a = a(1);
                    disp(sprintf('Using mode with neff=%.6f', m{a}.neff))
                    NewModes{mc} = m{a};
                    verbastelt = [ verbastelt a ];
                end

            else
                
                % Now find new modes based on evaluating the mode-field product E x H* and
                % finding the closest matches.
                                
                LastModes = Modes{end};
                for lmc = 1 : length(LastModes)
                    % This mode branch has been lost
                    if isempty(LastModes{lmc})
                        continue
                    end

                    % Find new mode which is closest to the old one
                    sprods = nan(1, length(m));
                    for i = 1 : length(m)
                        if length(find(verbastelt == i)) == 1
                            % New mode already used for a different branch
                            continue
                        end
                        
                        % Skalarprodukt von LastModes{lmc} mit m{i} bilden
                        % Alten Modus auf neues Grid interpolieren
                        [R,Z] = meshgrid(real(LastModes{lmc}.r), real(LastModes{lmc}.z));
                        laster = interp2(R, Z, LastModes{lmc}.er, real(m{i}.r), real(m{i}.z)', 'linear', 0.0);
                        lastez = interp2(R, Z, LastModes{lmc}.ez, real(m{i}.r), real(m{i}.z)', 'linear', 0.0);
                        a = (laster .* conj(m{i}.hz) - lastez .* conj(m{i}.hr)) / 2;
                        a = wgms3d_int(a, real(LastModes{lmc}.r), real(LastModes{lmc}.z)) ...
                            / (sqrt(m{i}.K) * sqrt(LastModes{lmc}.K));
                        sprods(i) = abs(a);
                    end
                    sprods
                    
                    % Den mit dem Maximum nehmen. Wenn Maximum < .950, dann noch mal mit
                    % kleinerer Stepsize. Wenn Stepsize zu klein wird,
                    % abbrechen.
                    maxsprods = max(sprods);
                    if maxsprods < 0.950
                        modelost = true;
                        stepok = 0;

                        assert(length(PARs) >= 1);
                        PAR = PARs(end);
                        PARstep = PARstep / 2;
                        
                        if abs(PARstep) >= abs(PARstep_min)
                            disp(sprintf('Lowering PARstep to %e and trying again', PARstep))
                            break
                        else
                            % Schrittweite zu klein geworden, aufhören.
                            disp(sprintf('PARstep became too low, stopping.'))
                            return
                        end
                    else
                        a = find(sprods == maxsprods); a = a(1);
                        disp(sprintf('Using mode with neff=%.6f and sprod=%.3f', m{a}.neff, maxsprods))
                        m{a}.sprod = maxsprods;
                        NewModes{lmc} = m{a};
                        verbastelt = [ verbastelt a ];
                    end
                end
            end

            if ~modelost
                ind = length(PARs)+1;
                PARs(ind) = PAR;

                for i = 1 : length(NewModes)
                    if ~isempty(NewModes{i})
                        Neffs(i, ind) = NewModes{i}.neff;
                    end
                end

                nsearch = abs(real(Neffs(:, end)));
                
                % Remember mode fields as well...
                Modes{ind} = NewModes;
                clear NewModes;
                clear LastModes;
                
                % If tracing went fine for some time, try to increase step size again.
                if stepok < 3
                    stepok = stepok + 1;
                else
                    if abs(PARstep) < abs(PARstep_orig)
                        PARstep = 2 * PARstep;
                        if abs(PARstep) > abs(PARstep_orig)
                            PARstep = PARstep_orig;
                        end
                        stepok = 2;
                    end
                end
            end

        else
            
            ind = length(PARs)+1;
            PARs(ind) = PAR;
            
        end

        if stepok > 0 & ~isempty(plotfunc)
            plotfunc(PARs, Modes);
            drawnow
        end

        wgms3d_tracemodes_update_parameter;
        if finished
            return
        end
            
    end

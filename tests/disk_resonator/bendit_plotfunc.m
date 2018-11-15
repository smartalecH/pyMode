
function bendit_plotfunc (Rcurvs, Modes)
    
    wgms3d_tracemodes_plotneffs (Rcurvs, Modes);
    
    for i = 1 : length(Modes{end})
        ccfig(10+i); clf
        wgms3d_plot_ht(Modes{end}{i}, 'Geometry', 'prkna.mgp', 'LogContours', 3);
        title(sprintf('Mode %d', i))
    end
    
function ccfig (h)
  
    if ishandle(h)
        set(0,'CurrentFigure',h)
    else
        figure(h)
    end

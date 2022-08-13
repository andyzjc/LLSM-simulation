function ClosePlots
    global Figures
    if isfield(Figures,'fig1') 
        delete(Figures.fig1)
    end
    if isfield(Figures,'fig2') 
        delete(Figures.fig2)
    end
    if isfield(Figures,'fig3') 
        delete(Figures.fig3)
    end
    
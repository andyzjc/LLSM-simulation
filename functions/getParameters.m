function getParameters(handles)
    global Data
    Data.NAmax = str2num(handles.edit_NAmax.String);
    Data.k_NAmax = Data.NAmax /Data.n * Data.k_wave_exc; % k
    Data.NAmin = str2num(handles.edit_NAmin.String);  
    Data.k_NAmin = Data.NAmin /Data.n * Data.k_wave_exc;
    Data.NAideal = (Data.NAmin + Data.NAmax)/2;
    Data.k_ideal = Data.k_wave_exc * Data.NAideal / Data.n;
    Data.dither_period = str2num(handles.edit_dither_period.String); % um
    Data.dither_step = str2num(handles.edit_dither_step.String); % number of s teps per dither period 
    Data.gauss_bound_width = str2num(handles.edit_gaussian_boundwidth.String); % Gaussian Bounding, um
    
    Data.angle = str2num(handles.edit_beam_angle.String);
    Data.weighting =  exp(str2num(handles.edit_beam_weighting.String));
    Data.apertureNA = str2num(handles.edit_gauss_circular_NA.String);
    Data.k_apertureNA = Data.apertureNA * Data.k_wave_exc / Data.n;
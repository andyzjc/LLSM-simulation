function setDefault(handles)
    global Data
    Data.sti_case = 'hex';

    handles.edit_NAmax.String = num2str(0.6);
    handles.edit_NAmin.String = num2str(0.5);
    handles.edit_dither_period.String = num2str(10);
    handles.edit_dither_step.String = num2str(201);
    handles.edit_gaussian_boundwidth.String = num2str(4);
    handles.edit_beam_angle.String = num2str([30, 90, 150, 210, 270, 330]);
    handles.edit_beam_weighting.String = num2str([0,0,0,0,0,0]);
    handles.edit_gauss_circular_NA.String = num2str(0.3);
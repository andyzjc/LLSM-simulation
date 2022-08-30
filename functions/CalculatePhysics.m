function Calculate_physics(handles)
    global Data

    % Physical Parameter 
    Data.N = str2num(handles.edit_Image_Size.String); % pixels
    Data.k_xz_scale = str2num(handles.edit_XZ_sampling.String);
    Data.y_scale = str2num(handles.edit_propagation_step.String);
    Data.n = str2num(handles.edit_refractive_index.String);
    Data.lambda_exc = str2num(handles.edit_wavelength_exc.String) / 1000; % um 
    Data.lambda_det = str2num(handles.edit_wavelength_det.String) / 1000;
    Data.wavelength_exc =  Data.lambda_exc / Data.n;
    Data.wavelength_det = Data.lambda_det / Data.n;

    [ax, ~] = meshgrid(  -(Data.N-1)/2 : (Data.N-1)/2 ) ; 
    Data.k_wave_exc = 1/Data.wavelength_exc;
    Data.k_bound = Data.k_xz_scale * Data.k_wave_exc; 
    Data.deltak = 2 * Data.k_bound / Data.N;
    Data.deltax = 1/(2 * Data.k_bound);

    % excitation
    % xz 
    Data.kx_exc = Data.deltak * ax;  %in unit wavelength
    Data.kz_exc = Data.kx_exc';
    Data.ky_exc = sqrt(Data.k_wave_exc^2 - Data.kx_exc.^2 - Data.kz_exc.^2);
    Data.ky_exc(Data.kx_exc.^2 + Data.kz_exc.^2 > Data.k_wave_exc.^2 ) = 0;
    Data.x_exc = Data.deltax * ax; 
    Data.z_exc = Data.x_exc'; 
    % y-propagation
    Data.y_exc = (  -(Data.N-1)/2 : (Data.N-1)/2  ) * Data.deltax * Data.y_scale; 
    Data.KY_exc = (-(Data.N-1)/2 : (Data.N-1)/2) * 1/(2*max(Data.y_exc)) / (2*Data.k_wave_exc);

    % for displaying
    Data.KX_exc = Data.kx_exc(1,:) / (2*Data.k_wave_exc);
    Data.KZ_exc = Data.KX_exc';
    Data.X_exc = Data.x_exc(1,:)  / Data.wavelength_exc; % value * wavelength = physical value (um)
    Data.Z_exc = Data.X_exc'; 
    Data.Y_exc = Data.y_exc / Data.wavelength_exc; 
 
    % detection
    Data.k_wave_det = 1/Data.wavelength_det;
    Data.NAdet = str2num(handles.edit_NAdet.String);
    Data.k_det = Data.k_wave_det * Data.NAdet / Data.n;
    %xy
    Data.kx_det = Data.deltak * ax;
    Data.ky_det = Data.kx_det';
    Data.kz_det = sqrt(Data.k_wave_det^2 - Data.kx_det.^2 - Data.ky_det.^2) .* (Data.k_det.^2 >= Data.kx_det.^2 + Data.ky_det.^2);
    Data.kz_det(Data.kx_det.^2 + Data.ky_det.^2 > Data.k_wave_det^2) = 0;
    Data.x_det = Data.deltax * ax;
    Data.y_det = Data.x_det';

    %z-propagation 
    Data.z_det = (  -(Data.N-1)/2 : (Data.N-1)/2  )  * Data.deltax * Data.y_scale; % no scaling factor
    Data.KZ_det = (-(Data.N-1)/2 : (Data.N-1)/2) * 1/(2*max(Data.z_det)) / (2*Data.k_wave_det);

    % for displaying
    Data.KX_det = Data.kx_det(1,:) / (2 * Data.k_wave_det);
    Data.KY_det = Data.KX_det';
    Data.X_det = Data.x_det(1,:)  / Data.wavelength_det;
    Data.Y_det = Data.X_det';
    Data.Z_det = Data.z_det / Data.wavelength_det;
    
    % Calculate FOV
    handles.text_step_size.String = num2str(Data.deltax,'%.2f') + " x " + ...
                                    num2str(Data.deltax,'%.2f') + " x " + ...
                                    num2str(Data.deltax * Data.y_scale,'%.2f');

    handles.text_FOV.String = num2str(Data.deltax*Data.N,'%.1f') + " x " + ...
                              num2str(Data.deltax* Data.N,'%.1f') + " x " + ...
                              num2str(max(Data.y_exc),'%.1f');

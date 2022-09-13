function Simulate(handles)
    global Data
    global Lattice
    tic
    % Define pupil function
    Data.Pupil_fun_exc = zeros(Data.N,Data.N);
    switch Data.sti_case
        case 'sw'
            CreateLattice;
            Lattice.islattice = 1;
        case 'square'
            CreateLattice;
            Lattice.islattice = 1;
        case 'hex'
            CreateLattice;
            Lattice.islattice = 1;
        case 'bessel'
            Data.Pupil_fun_exc = ((Data.k_NAmax > sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)) ...
                             .* (Data.k_NAmin < sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)));
            Lattice.islattice = 0;
        case 'other'
            CreateLattice;
            Lattice.islattice = 1;
        case '2dgaussian'
            Data.Pupil_fun_exc = exp( -(Data.kx_exc.^2 + Data.kz_exc.^2)/ ((Data.k_apertureNA/2).^2) );
            Lattice.islattice = 0;
        case '1dgaussian'
            Data.Pupil_fun_exc(:,(Data.N+1)/2) = exp( -(Data.kz_exc(:,1).^2)/ ((Data.k_apertureNA/2).^2) );
            Lattice.islattice = 0;
        case '2dairy'
            Data.Pupil_fun_exc = Data.k_apertureNA.^2 > (Data.kx_exc.^2 + Data.kz_exc.^2);
            Lattice.islattice = 0;
        case '1dairy'
            Data.Pupil_fun_exc(:,(Data.N+1)/2) = (Data.k_apertureNA) >= abs(Data.kz_exc(:,1));
            Lattice.islattice = 0;
    end

    disp("Case: " + Data.sti_case)

    % weighting to correct high NA 
    Data.Pupil_fun_det = Data.k_wave_exc./Data.kz_det;
    Data.Pupil_fun_det(Data.Pupil_fun_det == Inf) = 0;
    Data.Pupil_fun_det = fillmissing(Data.Pupil_fun_det,'constant',0);

    Data.Pupil_fun_exc = Data.Pupil_fun_exc .* Data.k_wave_exc./Data.ky_exc;
    Data.Pupil_fun_exc(Data.Pupil_fun_exc == inf) = 0;
    Data.Pupil_fun_exc = fillmissing(Data.Pupil_fun_exc,'constant',0);

    % Simulation start here
    PSF_exc_3d = zeros(Data.N,Data.N,Data.N);
    PSF_exc_3d_dither = PSF_exc_3d;
    PSF_det_3d = PSF_exc_3d;

    disp("Propagating")
    % propagation
    for i = 1:length(Data.y_exc)
        propagator_exc = exp(2*pi * 1i * Data.ky_exc * Data.y_exc(i));
        PSF_exc_3d(:,:,i) = abs( fftshift(ifft2(Data.Pupil_fun_exc .* propagator_exc) ) ).^2;

        propagator_det = exp(2*pi * 1i * Data.kz_det * Data.z_det(i));
        PSF_det_3d(:,:,i) = abs( fftshift( ifft2(Data.Pupil_fun_det .* propagator_det) ) ).^2;
    end  
    
    % dithering for lattice
    disp("Dithering")

    if handles.checkbox_fastdither.Value == 0
        if Lattice.islattice == 1 
            % dithering along x exc
            for j = 1:Data.dither_step
                PSF_exc_3d_dither = PSF_exc_3d_dither + ...
                    circshift(PSF_exc_3d,round(j * Data.dither_period / Data.deltax / Data.dither_step),2);
            end
        else
            for j = 1:length(Data.X_exc)*2
                PSF_exc_3d_dither = PSF_exc_3d_dither + ...
                    circshift(PSF_exc_3d,j,2);
            end
            Data.dither_period = max(Data.X_exc) *2;
            Data.dither_step = length(Data.X_exc) *2;
        end
    else
        for i = 1:size(PSF_exc_3d,3)
            PSF_exc_3d_dither(:,:,i) = meshgrid(mean(squeeze(PSF_exc_3d(:,:,i)),2))';
        end
        Data.dither_period = 0;
        Data.dither_step = 0;
    end

    % Normalize
    Data.PSF_exc_3d = PSF_exc_3d/max(max(max(PSF_exc_3d)));
    Data.PSF_exc_3d_dither = PSF_exc_3d_dither/max(max(max(PSF_exc_3d_dither)));
    Data.PSF_det_3d = PSF_det_3d/max(max(max(PSF_det_3d)));  
    toc
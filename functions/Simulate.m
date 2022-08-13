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
            Lattice.islattice = 1;
        case 'other'
            CreateLattice;
            Lattice.islattice = 1;
        case '2dgaussian'
            Data.Pupil_fun_exc = exp( -(Data.kx_exc.^2 + Data.kz_exc.^2)/ (Data.k_apertureNA) );
            Lattice.islattice = 0;
        case '1dgaussian'
            Data.Pupil_fun_exc(:,(Data.N+1)/2) = exp( (Data.KZ_exc.^2)/ (Data.k_apertureNA) );
            Lattice.islattice = 0;
        case '2dairy'
            Data.Pupil_fun_exc = Data.k_apertureNA.^2 > (Data.kx_exc.^2 + Data.kz_exc.^2);
            Lattice.islattice = 0;
        case '1dairy'
            Data.Pupil_fun_exc(:,(Data.N+1)/2) = Data.k_apertureNA >= (Data.KZ_exc.^2) ;
            Lattice.islattice = 0;
    end

    disp("Case: " + Data.sti_case)

    % Detection
    Data.Pupil_fun_det = Data.k_det.^2 > Data.kx_det.^2 + Data.ky_det.^2;

    % Simulation start here
    PSF_exc_3d = zeros(Data.N,Data.N,Data.N);
    PSF_exc_3d_dither = PSF_exc_3d;
    PSF_det_3d = PSF_exc_3d;

    disp("Propagating")
    % propagation
    for i = 1:length(Data.y_exc)
        propagator_exc = exp(2*pi * 1i * Data.ky_exc * Data.y_exc(i));
        PSF_exc_3d(:,:,i) = abs( fftshift( ifft2(Data.Pupil_fun_exc .* propagator_exc) ) ).^2;

        propagator_det = exp(2*pi * 1i * Data.kz_det * Data.z_det(i));
        PSF_det_3d(:,:,i) = abs( fftshift( ifft2(Data.Pupil_fun_det .* propagator_det) ) ).^2;
    end  
    
    % dithering for lattice
    disp("Dithering")

    % dithering along x exc
    for j = 1:Data.dither_step
        PSF_exc_3d_dither = PSF_exc_3d_dither + ...
            circshift(PSF_exc_3d,round(j * Data.dither_period / Data.deltax / Data.dither_step),2);
    end

    % Overall PSF
    Overall_PSF_axial = squeeze(PSF_exc_3d_dither(:,:,(Data.N+1)/2)) .* squeeze(PSF_det_3d(:,(Data.N+1)/2,:))' ;
    Overall_PSF_lateral = squeeze(PSF_exc_3d_dither(:,(Data.N+1)/2,:)) .* squeeze(PSF_det_3d(:,:,(Data.N+1)/2));
    
    % Normalize
    Data.PSF_exc_3d = PSF_exc_3d/max(max(max(PSF_exc_3d)));
    Data.PSF_exc_3d_dither = PSF_exc_3d_dither/max(max(max(PSF_exc_3d_dither)));
    Data.PSF_det_3d = PSF_det_3d/max(max(max(PSF_det_3d)));  
    Data.Overall_PSF_axial = Overall_PSF_axial/max(max(Overall_PSF_axial));
    Data.Overall_PSF_lateral = Overall_PSF_lateral/max(max(Overall_PSF_lateral));
    toc
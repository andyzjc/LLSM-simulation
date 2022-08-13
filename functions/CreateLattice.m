function CreateLattice
    global Lattice
    global Data
    N = Data.N;
    theta = Data.angle;
    weighting = Data.weighting;
    
    Lattice.Illumi_ideal = zeros(size(Data.x_exc));
    kxposition = Data.k_ideal * cosd(theta) /Data.deltak; % pixel
    kzposition = Data.k_ideal * sind(theta) /Data.deltak; % pixel
    
    for j = 1:length(kxposition)
    
        Lattice.Illumi_ideal( ...
            (N+1)/2 + round(kzposition(j)) ,...
            (N+1)/2 + round(kxposition(j)) ) = 1 * weighting(j);
    end
    Lattice.E_ideal = ifft2(Lattice.Illumi_ideal); 
    
    % bounded lattice 
    gauss_bound = exp(-2 * Data.z_exc.^2 / (Data.gauss_bound_width)^2);
    Lattice.E_bound = gauss_bound .* Lattice.E_ideal;
    
    % bounded back pupil
    Lattice.Illum_bound = abs(fft2(fftshift(Lattice.E_bound))).^2;
    Lattice.Illum_bound = Lattice.Illum_bound/max(max(Lattice.Illum_bound));
    
    % Generate mask 
    Lattice.A_mask = ((Data.k_NAmax > sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)) .* (Data.k_NAmin < sqrt(Data.kx_exc.^2 + Data.kz_exc.^2)));
    
    % Pupil functions
    Data.Pupil_fun_exc = Lattice.Illum_bound .* Lattice.A_mask;
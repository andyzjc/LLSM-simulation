clear all
close all

%% Desired Lattice
theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [30, 150, 210, 330]; % square
% theta = [90, 270]; % standing wave

%% Propagation in y-direction

% Physical Parameter 
Nxz = 1025; % pixels
Ny = 1025; % pixels
n = 1.33;
wavelength_exc = 0.488 / n; % um 
wavelength_dect = 0.488 / n;
NAdect = 1.1;
NAmin = 0.50;
NAmax = 0.60;
NAideal = (NAmin + NAmax)/2;
dither_period = 5; % um
dither_step = 201; % number of s teps per dither period 
gauss_bound_width = 1; % Gaussian Bounding, um
y_scale = 10; 

% when 2 plane wave traveling in opposite direction and interfere with each
% other, the periodicity of resulting standing wave is wavelength/4pi

% Real space image parameter
deltax = wavelength_exc/4/pi; 
deltaz = deltax;
deltay = deltax * y_scale;
deltax_dect = wavelength_dect / (2 * NAdect); % diffraction limit 
deltay_dect = deltax_dect;

% Fourier space image parameter
plot_scale = 1.5; % for ploting
K_exc = 1/(2*deltax); % Nyquist
K_dect = 1/(2*deltax_dect); 
deltakx = 2 * K_exc / Nxz; 
deltakz = deltakx;
deltakx_dect = 2 * K_dect / Nxz; 
deltakz_dect = deltakx_dect;
k_exc_ideal = NAideal/ (2 * pi) * K_exc;
k_dect_ideal = NAdect/ (2 * pi) * K_dect;

% Create mesh 
[ax, az] = meshgrid(  -(Nxz-1)/2 : (Nxz-1)/2 ) ; 
kx = deltakx * ax;  %in unit 1/um
kz = deltakz * az;
ky = sqrt(K_exc^2 - kx.^2 - kz.^2);
KX = kx(1,:) * wavelength_exc; % value / wavelength = k value (1/um)
KZ = kz(:,1)' * wavelength_exc;
x = deltax * ax;
z = deltaz * az; 
y = (0:Ny) * deltay; % max(y) = FOVy 
X = x(1,:) / wavelength_exc; % value * wavelength = physical value (um)
Z = z(:,1) / wavelength_exc; 
Y = y / wavelength_exc; 

% detection
x_dect = deltax_dect * ax;  
y_dect = x_dect';
Xdect = x_dect(1,:) / wavelength_exc; % pixels * wavelength = physical space
Ydect = y_dect(:,1) / wavelength_exc; 

kxdect = deltakx_dect * ax;
kydect = kxdect'; 
KXdect = kxdect(1,:) * wavelength_exc;
KYdect = kydect(:,1) * wavelength_exc;

%% Generate Ideal lattice 
Illumi_ideal = zeros(size(kx));

kxposition = k_exc_ideal * cosd(theta) / deltakx; % pixel
kzposition = k_exc_ideal * sind(theta) / deltakz; % pixel

% for ploting
KX_plot = floor((Nxz+1)/2 + min(kzposition)*plot_scale : ...
          (Nxz+1)/2 + max(kzposition)*plot_scale);
KZ_plot = KX_plot;
OTF_KX_plot = floor( (Nxz+1)* 6/16 : (Nxz+1) * 10/16 );
OTF_KZ_plot = OTF_KX_plot;

% Ideal lattice illumination
for j = 1:length(kxposition)

    Illumi_ideal( ...
        (Nxz+1)/2 + round(kzposition(j)) -1: (Nxz+1)/2 + round(kzposition(j)) +1,...
        (Nxz+1)/2 + round(kxposition(j)) ) = 1;
end

% Ideal lattice
E_ideal = fftshift(fft2(Illumi_ideal)); % no need to shift since already at center

fig1 = figure(1);
    subplot(3,4,1);

image1 = imagesc( KX(1,KX_plot),...
                  KZ(1,KZ_plot),...
                  Illumi_ideal(KX_plot,KZ_plot) );
    colormap(jet)
    title("Rear Pupil Illumination of ideal Lattice, " +...
          "N_{xz} = " + num2str(Nxz) + ...
          ", K_{bound} = " + num2str(K_exc))
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;

    subplot(3,4,2);
image2 = imagesc(X, Z ,abs(E_ideal));
    title("Ideal Lattice," + ...
          "\lambda_{exc}/n = " + num2str(wavelength_exc, '%.3f') + "um" + ...
          " , n = " + num2str(n))
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

%% Gaussian Bounding and get rear pupil illum back

gauss_bound = exp(-2 * z.^2 / (gauss_bound_width)^2);
E_bound = gauss_bound .* E_ideal;
    subplot(3,4,3)
image3 = imagesc(X, Z, abs(E_bound));
    title("Bounded Ideal Lattice, " + ...
          "Bound width = " + num2str(gauss_bound_width) + "um")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

Illum_bound = abs(ifft2(ifftshift(E_bound))).^2;
Illum_bound = Illum_bound/max(max(Illum_bound));
    subplot(3,4,4)
image4 = imagesc( KX(1,KX_plot),...
                  KZ(1,KZ_plot),...
                  Illum_bound(KX_plot,KZ_plot) );
    title("Rear pupil illumination after bounding")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;

%% Generate mask and filter
k_NAmax = NAmax/wavelength_exc; % k
k_NAmin = NAmin/wavelength_exc; 

% create mask
A_mask = ((k_NAmax > sqrt(kx.^2 + kz.^2)) .* (k_NAmin < sqrt(kx.^2 + kz.^2)));
Illum_mask = imfuse(Illum_bound,A_mask,"falsecolor","ColorChannels","green-magenta");
    subplot(3,4,5)
image5 = imagesc( KX(1,KX_plot),...
                  KZ(1,KZ_plot),...
                  Illum_mask(KX_plot,KZ_plot,:) );
    title("Masking, " +...
          "NA_{max} = " + num2str(NAmax) +...
          ", NA_{min} = " + num2str(NAmin) )
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image;

% Create pupil function    
Pupil_fun_exc = Illum_bound .* A_mask;
    subplot(3,4,6)
image6 = imagesc( KX(1,KX_plot),...
                  KZ(1,KZ_plot),...
                  Pupil_fun_exc(KX_plot,KZ_plot) );
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;


%% PSF and OTF
PSF_exc = abs( fftshift(fft2(Pupil_fun_exc)) ).^2;
PSF_exc = PSF_exc/max(max(PSF_exc));
    subplot(3,4,7)
image7 = imagesc(X, Z, PSF_exc);
    title("Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

% OTF
OTF_exc = abs(fftshift(fft2(PSF_exc)));
OTF_exc = OTF_exc/max(max(OTF_exc));
    subplot(3,4,8)
image8 = imagesc( KX( OTF_KX_plot ),...
                  KZ( OTF_KZ_plot ),...
                  OTF_exc( OTF_KX_plot, OTF_KZ_plot ) ) ;
    title("Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;

%% Dither 
PSF_exc_dither = PSF_exc;
for j = 1:dither_step
    PSF_exc_dither = PSF_exc_dither + ...
                     circshift(PSF_exc,round(j * dither_period / deltax / dither_step),2)/2;
end
PSF_exc_dither = PSF_exc_dither / max(max(PSF_exc_dither));

    subplot(3,4,9)
    hold on;
image9 = imagesc(X, Z ,PSF_exc_dither);
    title("Dithered Excitation PSF, " + ...
          "T_d = " + num2str(dither_period) + "um, " +...
          "Dither Step = " + num2str(dither_step))
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;
line1 = xline(X((Nxz+1)/2));    
    line1.Color = 'r';
    line1.LineWidth = 1;
    line1.LineStyle = '--';
    hold off

    subplot(3,4,10)
PSF_exc_dither_profile = PSF_exc_dither(:,(Nxz+1)/2);
PSF_exc_dither_profile = PSF_exc_dither_profile/max(max(PSF_exc_dither_profile));
image10 = plot( Z(OTF_KZ_plot), PSF_exc_dither_profile(OTF_KZ_plot));
    title("PSF Dithered Line Profile")
    ylabel("Normalized a.u. ")
    xlabel("z/\lambda")
    image10.LineWidth = 2;
    image10.Color = 'r';
    image10.Parent.YAxis.TickValues = linspace(0,1,11);
    grid on

%% Dither OTF 
OTF_dither = abs(fftshift(fft2(PSF_exc_dither)));
OTF_dither = OTF_dither/max(max(OTF_dither));

    subplot(3,4,11)
    hold on
image11 = imagesc( KX(OTF_KX_plot),...
                   KZ(OTF_KZ_plot),...
                   OTF_dither( OTF_KX_plot,...
                               OTF_KX_plot )  );
    title("Dithered Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image; 

    subplot(3,4,12)
OTF_dither = OTF_dither(:,(Nxz+1)/2);
OTF_dither = OTF_dither/max(max(OTF_dither));
image12 = plot( KZ(OTF_KZ_plot), OTF_dither(OTF_KZ_plot));
    title("OTF Dithered Line Profile, " + ...
        "K_x = " + num2str(0) )
    ylabel("Normalized a.u. ")
    xlabel("kz * \lambda")
    image12.LineWidth = 2;
    image12.Color = 'r';
    image12.Parent.YAxis.TickValues = linspace(0,1,11);
    grid on

%% Propagator
tic
E_prop = zeros(Nxz,Nxz,Ny);
I_prop = E_prop;
Ein = Pupil_fun_exc;
for i = 1:length(y)
    propagator = exp(1i * ky * y(i));
    Eout = fftshift(fft2((Ein .* propagator)));
    E_prop(:,:,i) = Eout;
    I_prop(:,:,i) = abs(Eout.^2);
end  
I_prop = I_prop/max(max(max(I_prop)));
toc

%% Detection
Pupil_fun_dect = k_dect_ideal > sqrt(kxdect.^2 + kydect.^2);

% detection PSF
PSF_dect = abs( fftshift(fft2(Pupil_fun_dect)) ).^2;
PSF_dect = PSF_dect/max(max(PSF_dect));
PSF_dect_profile = PSF_dect(:,(Nxz+1)/2);
PSF_dect_profile = PSF_dect_profile/max(max(PSF_dect_profile));

% detection OTF
OTF_dect = abs(fftshift(fft2(PSF_dect)));
OTF_dect = OTF_dect/max(max(OTF_dect));
OTF_dect_profile = OTF_dect(:,(Nxz+1)/2);
OTF_dect_profile = OTF_dect_profile/max(max(OTF_dect_profile));

%% Ploting
fig2 = figure(2);
    colormap(jet)

    % plot xz
    subplot(3,4,1:2);
    slice = 1; % focal plane
image13 = imagesc(X,Z, squeeze(I_prop(:,:,slice)));
    title("xz plane, " + "Y = " + num2str(y(slice)) +...
          ", NA_{max} = " + num2str(NAmax) + ...
          ", NA_{min} = " + num2str(NAmin) )
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar

    % plot yz profile
    subplot(3,4,3)
    slice = (Nxz+1)/2;
image14 = plot(Y,squeeze(I_prop(slice,slice,:)));
    title("PSF-yz, " + ...
          "X = " + num2str(X(slice)) + ...
          ", Z = " + num2str(Z(slice))  )
    xlabel("y/\lambda")
    ylabel("Normalized a.u. ")
    image14.Color = 'r';
    image14.LineWidth = 2;
    grid on

    % plot detection pupil func
    subplot(3,4,4)
image15 = imagesc(KXdect, KYdect, Pupil_fun_dect);
    title("Detection PSF, " + "NA = " + num2str(NAdect))
    xlabel("kx * \lambda")
    ylabel("ky * \lambda")
    colorbar
    axis image;

    % plot yz
    subplot(3,4,5:6)
    slice = (Nxz+1)/2; 
image16 = imagesc(Y,Z, squeeze(I_prop(:,slice,:)));
    title("yz plane, " + "X = " + num2str(X(slice)) )
    xlabel("y/\lambda")
    ylabel("z/\lambda")
    colorbar

    % plot detection PSF
    subplot(3,4,7)
image17 = imagesc(Xdect, Ydect, PSF_dect);
    title("Detection PSF")
    xlabel("x/\lambda")
    ylabel("y/\lambda")
    colorbar;
    axis image;  
    image17.Parent.XLim = [-10,10];
    image17.Parent.YLim = [-10,10];

    % plot PSF line profile
    subplot(3,4,8)
image18 = plot(Ydect, PSF_dect_profile);
    title("PSF_{dect} profile, " + ...
          "X = " + num2str(0) )
    xlabel("y/\lambda")
    ylabel("Normalized a.u. ")
    image18.Color = 'r';
    image18.LineWidth = 2;
    image18.Parent.XLim = [-10,10];
    grid on
    
    % plot xy
    subplot(3,4,9:10)
    slice = (Nxz+1)/2;
image19 = imagesc(Y,X,squeeze(I_prop(slice,:,:)));
    title("xy plane, " + "Z = " + num2str(Z(slice)))
    xlabel("y/\lambda")
    ylabel("x/\lambda")
    colorbar

    % plot detection OTF
    subplot(3,4,11)
image17 = imagesc(KXdect, KYdect, OTF_dect);
    title("Detection OTF")
    xlabel("kx * \lambda")
    ylabel("ky * \lambda")
    colorbar;
    axis image;  

    % plot OTF line profile
    subplot(3,4,12)
image18 = plot(KYdect, OTF_dect_profile);
    title("OTF_{dect} profile, " + ...
          "KX = " + num2str(0) )
    xlabel("ky * \lambda")
    ylabel("Normalized a.u. ")
    image18.Color = 'r';
    image18.LineWidth = 2;
    image18.Parent.XLim = [-0.5 0.5];
    grid on



    


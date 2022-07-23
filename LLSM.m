close all
clear all

%% Desired Lattice
theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [30, 150, 210, 330]; % square
% theta = [90, 270]; % standing wave

% Physical Parameter 
Nxz = 1025; % pixels
Ny = 257; % pixels
n = 1.33;
wavelength = 0.488 / n; % um 
NAdect = 1;
NAmin = 0.57;
NAmax = 0.65;
NAideal = (NAmin + NAmax)/2;
dither_period = 5; % um
dither_step = 201; % number of steps per dither period 
a = 5 ; % Gaussian Bounding, um

% when 2 plane wave traveling in opposite direction and interfere with each
% other, the periodicity of resulting standing wave is wavelength/4pi

% Real space image parameter
deltax = wavelength/4/pi; 
deltaz = deltax;

% Fourier space image parameter
plot_scale = 1.5; % for ploting
Kbound = 1/(2*deltax); % Nyquist
deltakx = 2*Kbound / Nxz; 
deltakz = deltakx;
k_ideal = NAideal/ ( 2* pi) * Kbound;

% Create mesh 
[ax, az] = meshgrid(  -(Nxz-1)/2 : (Nxz-1)/2 ) ; 
kx = deltakx * ax;  %in unit 1/um
kz = deltakz * az;
x = 4 * pi * deltax * ax;
z = 4 * pi * deltaz * az; 
X = x(1,:); % pixel length = wavelength
Z = z(:,1); 
KX = kx(1,:); % pixel length = 1/wavelength
KZ = kz(:,1)';

%% Generate Ideal lattice 
Illumi_ideal = zeros(size(kx));

kxposition = k_ideal * cosd(theta) / deltakx; % pixel
kzposition = k_ideal * sind(theta) / deltakz; % pixel

% for ploting
KX_plot = floor((Nxz+1)/2 + min(kzposition)*plot_scale : ...
          (Nxz+1)/2 + max(kzposition)*plot_scale);
KZ_plot = KX_plot;
OTF_KX_plot = floor( (Nxz+1)* 6/16 : (Nxz+1) * 10/16 );
OTF_KZ_plot = OTF_KX_plot;


% Ideal lattice illumination
for j = 1:length(kxposition)

    Illumi_ideal( ...
        (Nxz+1)/2 + round(kzposition(j)) ,...
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
          ", K_{bound} = " + num2str(Kbound))
    xlabel("kx(4\pi/\lambdaNxz) (1/um)")
    ylabel("kz(4\pi/\lambdaNxz) (1/um)")
    colorbar;
    axis image;

    subplot(3,4,2);
image2 = imagesc(X, Z ,abs(E_ideal));
    title("Ideal Lattice, " + ...
          "\lambda_{exc}/n = " + num2str(wavelength, '%.3f') + "um" + ...
          " , n = " + num2str(n))
    xlabel("x(\lambda) (um)")
    ylabel("z(\lambda) (um)")
    colorbar;
    axis image;

%% Gaussian Bounding and get rear pupil illum back

gauss_bound = exp(-2 * z.^2 / a^2);
E_bound = gauss_bound .* E_ideal;
    subplot(3,4,3)
image3 = imagesc(X, Z, abs(E_bound));
    title("Bounded Ideal Lattice, " + ...
          "Bound width = " + num2str(a) + "um")
    xlabel("x(\lambda) (um)")
    ylabel("z(\lambda) (um)")
    colorbar;
    axis image;

Illum_bound = abs(ifft2(ifftshift(E_bound))).^2;
Illum_bound = Illum_bound/max(max(Illum_bound));
    subplot(3,4,4)
image4 = imagesc( KX(1,KX_plot),...
                  KZ(1,KZ_plot),...
                  Illum_bound(KX_plot,KZ_plot) );
    title("Rear pupil illumination after bounding")
    xlabel("kx(4\pi/\lambdaNxz) (1/um)")
    ylabel("kz(4\pi/\lambdaNxz) (1/um)")
    colorbar;
    axis image;

%% Generate mask and filter
k_NAmax = NAmax/wavelength; % k
k_NAmin = NAmin/wavelength; 

% create mask
A_mask = ((k_NAmax > sqrt(kx.^2 + kz.^2)) .* (k_NAmin < sqrt(kx.^2 + kz.^2)));
Illum_mask = imfuse(Illum_bound,A_mask);
    subplot(3,4,5)
image5 = imagesc( KX(1,KX_plot),...
                  KZ(1,KZ_plot),...
                  Illum_mask(KX_plot,KZ_plot) );
    title("Masking, " +...
          "NA_{max} = " + num2str(NAmax) +...
          ", NA_{min} = " + num2str(NAmin) )
    xlabel("kx(4\pi/\lambdaNxz) (1/um)")
    ylabel("kz(4\pi/\lambdaNxz) (1/um)")
    axis image;

% Create pupil function    
Pupil_fun = Illum_bound .* A_mask;
    subplot(3,4,6)
image6 = imagesc( KX(1,KX_plot),...
                  KZ(1,KZ_plot),...
                  Pupil_fun(KX_plot,KZ_plot) );
    title("Pupil Function")
    xlabel("kx(4\pi/\lambdaNxz) (1/um)")
    ylabel("kz(4\pi/\lambdaNxz) (1/um)")
    colorbar;
    axis image;


%% PSF and OTF
PSF_exc = abs( fftshift(fft2(Pupil_fun)) ).^2;
PSF_exc = PSF_exc/max(max(PSF_exc));
    subplot(3,4,7)
image7 = imagesc(X, Z, PSF_exc);
    title("Excitation PSF")
    xlabel("x(\lambda) (um)")
    ylabel("z(\lambda) (um)")
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
    xlabel("kx(4\pi/\lambdaNxz) (1/um)")
    ylabel("kz(4\pi/\lambdaNxz) (1/um)")
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
    xlabel("x(\lambda) (um)")
    ylabel("z(\lambda) (um)")
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
    xlabel("z(\lambda) (um)")
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
    xlabel("kx(4\pi/\lambdaNxz) (1/um)")
    ylabel("kz(4\pi/\lambdaNxz) (1/um)")
    colorbar;
    axis image; 

    subplot(3,4,12)
OTF_dither = OTF_dither(:,(Nxz+1)/2);
OTF_dither = OTF_dither/max(max(OTF_dither));
image12 = plot( Z(OTF_KZ_plot), OTF_dither(OTF_KZ_plot));
    title("OTF Dithered Line Profile")
    ylabel("Normalized a.u. ")
    xlabel("z(\lambda) (um)")
    image12.LineWidth = 2;
    image12.Color = 'r';
    image12.Parent.YAxis.TickValues = linspace(0,1,11);
    grid on

%% Propagation in y-direction
y_scale = 100; 
deltay = wavelength/4/pi * y_scale;
ky = sqrt(Kbound^2 - kx.^2 - kz.^2);
y = (0:Ny) * deltay; % max(y) = FOVy 

E_prop = zeros(Nxz,Nxz,Ny);
I_prop = E_prop;
Ein = Pupil_fun/(max(max(Pupil_fun)));

%% Propagator

% for i = 1:length(y)
%     propagator = exp(1i * ky * y(i));
%     Eout = fftshift(fft2((Ein .* propagator)));
%     imagesc(x(1,:)/wavelength,z(:,1)/wavelength,abs(Eout.^2))
%     drawnow
%     colormap(jet)
%     pause(0.5)
%     caxis([0 4000])
%     disp(i)
%     E_prop(:,:,i) = Eout;
%     I_prop(:,:,i) = abs(Eout.^2);
% end  




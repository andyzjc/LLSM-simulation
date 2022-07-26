clear all
close all

%% Desired Lattice
theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [30, 150, 210, 330]; % square
% theta = [90, 270]; % standing wave

% Physical Parameter 
Nxz = 513; % pixels
Ny = 513; % pixels
n = 1.33;
lambda_exc = 0.488; % um 
lambda_dect = 0.488;
wavelength_exc =  lambda_exc / n;
wavelength_dect = lambda_dect / n;
NAmin = 0.50;
NAmax = 0.60;
NAdect = 1.1;
NAideal = (NAmin + NAmax)/2;
dither_period = 5; % um
dither_step = 201; % number of s teps per dither period 
gauss_bound_width = 3; % Gaussian Bounding, um
xz_scale = 2;
y_scale = 1;

k_wave = 1/wavelength_exc;
k_ideal = k_wave * NAideal / n;
k_bound = xz_scale * k_wave;
deltak = 2 * k_bound / Nxz;
deltax = 1/(2 * k_bound);

[ax, az] = meshgrid(  -(Nxz-1)/2 : (Nxz-1)/2 ) ; 
kx = deltak * ax;  %in unit wavelength
kz = kx';
ky = sqrt(k_wave^2 - kx.^2 - kz.^2);
x = deltax * ax; 
z = x'; 
y = (0:Ny) * deltax * y_scale; % max(y) = FOVy 

KX = kx(1,:) / k_wave;
KZ = KX';

X = x(1,:)  / wavelength_exc; % value * wavelength = physical value (um)
Z = X'; 
Y = y / wavelength_exc; 

% Generate Ideal lattice 
Illumi_ideal = zeros(size(ax));

kxposition = k_ideal * cosd(theta) /deltak; % pixel
kzposition = k_ideal * sind(theta) /deltak; % pixel

% Ideal lattice illumination
for j = 1:length(kxposition)

    Illumi_ideal( ...
        (Nxz+1)/2 + round(kzposition(j)) ,...
        (Nxz+1)/2 + round(kxposition(j)) ) = 1;
end

% Ideal lattice
E_ideal = ifft2(fftshift(Illumi_ideal)); % no need to shift since already at center


     subplot(3,4,1);
image1 = imagesc(KX,KZ, Illumi_ideal );
    colormap(jet)
    title("Rear Pupil Illumination of ideal Lattice, " +...
          "N_{xz} = " + num2str(Nxz) + ...
          ", K_{bound} = " + num2str(k_bound))
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
%     image1.Parent.XLim = [-xz_scale/2, xz_scale/2];
%     image1.Parent.YLim = [-xz_scale/2, xz_scale/2];
    colorbar;    

    subplot(3,4,2);
image2 = imagesc(X,Z, abs( E_ideal));
    title("Ideal Lattice," + ...
          "\lambda_{exc}/n = " + num2str(wavelength_exc, '%.3f') + "um" + ...
          " , n = " + num2str(n))
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    

% Gaussian Bounding and get rear pupil illum back

gauss_bound = exp(-2 * z.^2 / (gauss_bound_width)^2);
E_bound = gauss_bound .* E_ideal;

    subplot(3,4,3);
image3 = imagesc(X, Z, abs(E_bound));
    title("Bounded Ideal Lattice, " + ...
          "Bound width = " + num2str(gauss_bound_width) + "um")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    

Illum_bound = abs(ifftshift(fft2(ifftshift(E_bound)))).^2;
Illum_bound = Illum_bound/max(max(Illum_bound));

    subplot(3,4,4)
image4 = imagesc( KX, KZ,...
                  Illum_bound );
    title("Rear pupil illumination after bounding")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
%     image4.Parent.XLim = [-xz_scale/2, xz_scale/2];
%     image4.Parent.YLim = [-xz_scale/2, xz_scale/2];
    colorbar;
    

% Generate mask and filter
k_NAmax = NAmax /n * k_wave; % k
k_NAmin = NAmin /n * k_wave; 

% create mask
A_mask = ((k_NAmax > sqrt(kx.^2 + kz.^2)) .* (k_NAmin < sqrt(kx.^2 + kz.^2)));
Illum_mask = imfuse(Illum_bound,A_mask,"falsecolor","ColorChannels","green-magenta");
    subplot(3,4,5)
image5 = imagesc( KX, KZ,...
                  Illum_mask);
    title("Masking, " +...
          "NA_{max} = " + num2str(NAmax) +...
          ", NA_{min} = " + num2str(NAmin) )
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image;

% Create pupil function    
Pupil_fun_exc = Illum_bound;
    subplot(3,4,6)
image6 = imagesc( KX, KZ,...
                  Pupil_fun_exc );
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;


% PSF and OTF
PSF_exc = abs( fftshift(ifft2(Pupil_fun_exc)) ).^2;
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
image8 = imagesc( KX,...
                  KZ,...
                  OTF_exc ) ;
    title("Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;

% Dither 
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

    subplot(3,4,10)
PSF_exc_dither_profile = PSF_exc_dither(:,(Nxz+1)/2);
PSF_exc_dither_profile = PSF_exc_dither_profile/max(max(PSF_exc_dither_profile));
image10 = plot( Z, PSF_exc_dither_profile);
    title("PSF Dithered Line Profile")
    ylabel("Normalized a.u. ")
    xlabel("z/\lambda")
    image10.LineWidth = 2;
    image10.Color = 'r';
    image10.Parent.YAxis.TickValues = linspace(0,1,11);
    grid on

% Dither OTF 
OTF_dither = abs(fftshift(fft2(PSF_exc_dither)));
OTF_dither = OTF_dither/max(max(OTF_dither));

    subplot(3,4,11)
    hold on
image11 = imagesc( KX,...
                   KZ,...
                   OTF_dither  );
    title("Dithered Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image; 

    subplot(3,4,12)
OTF_dither = OTF_dither(:,(Nxz+1)/2);
OTF_dither = OTF_dither/max(max(OTF_dither));
image12 = plot( KZ, OTF_dither);
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

I_prop = zeros(Nxz,Nxz,Ny);
Ein = Pupil_fun_exc;
for i = 1:length(y)
    propagator = exp(2*pi * 1i * ky * y(i));
    Eout = fftshift(fft2((Ein .* propagator)));
    I_prop(:,:,i) = abs(Eout.^2);
end  
I_prop = I_prop/max(max(max(I_prop)));
toc

% %% Detection
% Pupil_fun_dect = k_dect_ideal > sqrt(kxdect.^2 + kydect.^2);
% 
% % detection PSF
% PSF_dect = abs( fftshift(fft2(Pupil_fun_dect)) ).^2;
% PSF_dect = PSF_dect/max(max(PSF_dect));
% PSF_dect_profile = PSF_dect(:,(Nxz+1)/2);
% PSF_dect_profile = PSF_dect_profile/max(max(PSF_dect_profile));
% 
% % detection OTF
% OTF_dect = abs(fftshift(fft2(PSF_dect)));
% OTF_dect = OTF_dect/max(max(OTF_dect));
% OTF_dect_profile = OTF_dect(:,(Nxz+1)/2);
% OTF_dect_profile = OTF_dect_profile/max(max(OTF_dect_profile));

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
image14 = plot(y,squeeze(I_prop(slice,slice,:)));
    title("PSF-yz, " + ...
          "X = " + num2str(X(slice)) + ...
          ", Z = " + num2str(Z(slice))  )
    xlabel("y/\lambda")
    ylabel("Normalized a.u. ")
    image14.Color = 'r';
    image14.LineWidth = 2;
    grid on

%     % plot detection pupil func
%     subplot(3,4,4)
% image15 = imagesc(KXdect, Kbound, Pupil_fun_dect);
%     title("Detection PSF, " + "NA = " + num2str(NAdect))
%     xlabel("kx * \lambda")
%     ylabel("ky * \lambda")
%     colorbar
%     axis image;

    % plot yz
    subplot(3,4,5:6)
    slice = (Nxz+1)/2; 
image16 = imagesc(y,Z, squeeze(I_prop(:,slice,:)));
    title("yz plane, " + "X = " + num2str(X(slice)) )
    xlabel("y/\lambda")
    ylabel("z/\lambda")
    colorbar
    axis image

%     % plot detection PSF
%     subplot(3,4,7)
% image17 = imagesc(Xdect, Ydect, PSF_dect);
%     title("Detection PSF")
%     xlabel("x/\lambda")
%     ylabel("y/\lambda")
%     colorbar;
%     axis image;  
%     image17.Parent.XLim = [-10,10];
%     image17.Parent.YLim = [-10,10];
% 
%     % plot PSF line profile
%     subplot(3,4,8)
% image18 = plot(Ydect, PSF_dect_profile);
%     title("PSF_{dect} profile, " + ...
%           "X = " + num2str(0) )
%     xlabel("y/\lambda")
%     ylabel("Normalized a.u. ")
%     image18.Color = 'r';
%     image18.LineWidth = 2;
%     image18.Parent.XLim = [-10,10];
%     grid on
%     
%     % plot xy
%     subplot(3,4,9:10)
%     slice = (Nxz+1)/2;
% image19 = imagesc(Y,X,squeeze(I_prop(slice,:,:)));
%     title("xy plane, " + "Z = " + num2str(Z(slice)))
%     xlabel("y/\lambda")
%     ylabel("x/\lambda")
%     colorbar
% 
%     % plot detection OTF
%     subplot(3,4,11)
% image17 = imagesc(KXdect, KYdect, OTF_dect);
%     title("Detection OTF")
%     xlabel("kx * \lambda")
%     ylabel("ky * \lambda")
%     colorbar;
%     axis image;  
% 
%     % plot OTF line profile
%     subplot(3,4,12)
% image18 = plot(KYdect, OTF_dect_profile);
%     title("OTF_{dect} profile, " + ...
%           "KX = " + num2str(0) )
%     xlabel("ky * \lambda")
%     ylabel("Normalized a.u. ")
%     image18.Color = 'r';
%     image18.LineWidth = 2;
%     image18.Parent.XLim = [-0.5 0.5];
%     grid on

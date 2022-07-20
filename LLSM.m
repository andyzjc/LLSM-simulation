close all
clear all

%% Initialize some parameters for back pupil plane
Nxz = 1025; % pixels
Ny = 1025; % pixels

% Physical Parameter 
wavelength = 0.488; % um
n = 1.33; % water
NAdect = 1.0;
NAmin = 0.57;
NAmax = 0.65;
NAideal = (NAmin + NAmax)/2;

%NAideal =0.6;

% Real space image parameter
xz_scale = 2;
deltax = wavelength / (xz_scale * NAdect); 
deltaz = deltax;

% Fourier space image parameter
k0 = 2*pi/wavelength;
k0_med = k0 * 1.33;
deltakx = 2*pi / (Nxz * deltax);
deltakz = deltakx;
Kbound = Nxz/2 * deltakx;

% Create mesh 
[ax, az] = meshgrid(  -(Nxz-1)/2 : (Nxz-1)/2 ) ; 
kx = deltakx * ax;  %in unit 1/um
kz = deltakz * az;
x = deltax * ax;
z = deltaz * az; 

%% Generate Ideal lattice 
Illumi_ideal = zeros(size(kx));
theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [30, 150, 210, 330]; % square
theta = [90, 270]; % standing wave
k_ideal = NAideal / wavelength; % in unit of 1/um
kxposition = k_ideal * cosd(theta);
kzposition = k_ideal * sind(theta);

% Ideal lattice illumination
for j = 1:length(kxposition)
    Illumi_ideal(...
        ((Nxz-1)/2 + round(kzposition(j)/deltakz)):((Nxz-1)/2 + round(kzposition(j)/deltakz))+1, ...
        (Nxz-1)/2 + round(kxposition(j)/deltakx)) = 1;
end

% Ideal lattice
E_ideal = fftshift(fft2(Illumi_ideal));

fig1 = figure(1);
    image1 = imagesc(x(1,:)/wavelength, z(:,1)/wavelength,abs(E_ideal));
    colormap("jet")
    title("Ideal Lattice")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

%% Gaussian Bounding and get rear pupil illum back
a = 30*wavelength;
gauss_bound = exp(-2 * z.^2 / a^2);
E_bound = gauss_bound .* E_ideal;
fig2 = figure(2);
    image2 = imagesc(x(1,:)/wavelength, z(:,1)/wavelength,abs(E_bound));
    colormap("jet")
    title("Bounded Ideal Lattice")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

%% Reverse
Illum_bound = abs(fft2(fftshift(E_bound))).^2;
fig3 = figure(3);
    image3 = imagesc(kx(1,:)/Kbound, kz(:,1)/Kbound,Illum_bound);
    colormap(jet)
    title("Bounded Illumination")
    xlabel("kx (Normalized by 2*pi/lambda)")
    ylabel("kz (Normalized by 2*pi/lambda)")
    colorbar;
    axis image;

%% Generate mask and filter
A_mask = zeros(size(kx));
k_NAmax = NAmax/wavelength; % 1/um
k_NAmin = NAmin/wavelength; 
k_diff = k_NAmax - k_NAmin; 

% radius = k_ideal; % pixels
% thickness = k_diff; % pixels

A_mask = ((k_NAmax >= sqrt(kx.^2 + kz.^2)).* (k_NAmin <= sqrt(kx.^2 + kz.^2)));
fig4 = figure(4);
    image4 = imagesc(kx(1,:)/Kbound, kz(:,1)/Kbound, A_mask);
    title("Mask")
    xlabel("kx (Normalized by 2*pi/lambda)")
    ylabel("kz (Normalized by 2*pi/lambda)")
    colorbar;
    axis image;
Pupil_fun = Illum_bound .* A_mask;
fig5 = figure(5);
    image1 = imagesc(kx(1,:)/Kbound, kz(:,1)/Kbound, Pupil_fun);
    colormap(jet)
    title("Pupil Function")
    xlabel("kx (Normalized by 2*pi/lambda)")
    ylabel("kz (Normalized by 2*pi/lambda)")
    colorbar;
    axis image;


%% PSF(x,z)
PSF_exc = abs( fftshift(fft2(fftshift(Pupil_fun) )) ).^2;
fig5 = figure(6);
    image2 = imagesc(x(1,:)/wavelength, z(:,1)/wavelength,PSF_exc);
    caxis([min(min(PSF_exc)) max(max(PSF_exc))]);
    colormap(jet)
    title("Excitation PSF")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

%% OTF
OTF_exc = abs(fftshift(fft2(PSF_exc)));
fig6 = figure(6);
    image3 = imagesc(kx(1,:)/Kbound, kz(:,1)/Kbound,OTF_exc);
    caxis([min(min(OTF_exc)) max(max(OTF_exc))]);
    colormap(jet)
    title("Excitation OTF")
    xlabel("kx (Normalized by 2*pi/lambda)")
    ylabel("kz (Normalized by 2*pi/lambda)")
    colorbar;
    axis image;

%% Dither 
PSF_exc_dither = meshgrid(sum(PSF_exc,2))';  
plot(z(:,1),PSF_exc_dither(:,513)/max(PSF_exc_dither(:,513)))
fig7 = figure(7);
    image4 = imagesc(x(1,:)/wavelength, z(:,1)/wavelength,PSF_exc_dither);
    caxis([min(min(PSF_exc_dither)) max(max(PSF_exc_dither))]);
    colormap(jet)
    title("Excitation PSF - Dither")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

%% Dither OTF 
OTF_dither = abs(fftshift(fft2(PSF_exc_dither)));

fig8 = figure(8);
    image5 = imagesc(kx(1,:)/Kbound, kz(:,1)/Kbound, OTF_dither);
    caxis([min(min(OTF_dither)) max(max(OTF_dither))]);
    colormap(jet)
    title("Excitation OTF - Dither")
    xlabel("kx (Normalized by 2*pi/lambda)")
    ylabel("kz (Normalized by 2*pi/lambda)")
    colorbar;
    axis image;

%% Propagation in y-direction
y_scale = 1;
deltay = wavelength / n *y_scale; % 1 pixel = 50 wavelength
ky = sqrt(k0_med^2 - kx.^2 - kz.^2);
y = (0:Ny) * deltay; % max(y) = FOVy 

E_prop = zeros(Nxz,Nxz,Ny);
I_prop = E_prop;
Ein = Pupil_fun/(max(max(Pupil_fun)));

%% Propagator

transverse = exp(1i * kx * x + 1i * kz * z);
% for i = 1:length(y)
for i = 1:10:500
    propagator = exp(1i * ky * y(i));
    Eout = fftshift(fft2((Ein .* propagator)));
    imagesc(x(1,:)/wavelength,z(:,1)/wavelength,abs(Eout.^2))
    drawnow
    colormap(jet)
    pause(0.5)
    caxis([0 4000])
    E_prop(:,:,i) = Eout;
    I_prop(:,:,i) = abs(Eout.^2);
end  

%% Fresnel Propagation
% Ein = fftshift(fft2(Pupil_fun ));
% for i = 1:length(y)
%     h = 1/(1i * wavelength * y(i)) * exp(1i * k0_med/( 2* y(i)) * (x.^2+z.^2));
%     Eout = propIR(Ein,h);
%     imagesc(abs(Eout.^2))
%     drawnow
%     disp(i)
%     pause(0.25)
%     E_prop(:,:,i) = Eout;
%     I_prop(:,:,i) = abs(Eout.^2);
% end
% 
% function[Eout]=propIR(Ein,h)
% 
%     H = fft2(fftshift(h)); %create trans func
%     F_Ein = fft2(fftshift(Ein)); %shift, fft src field
%     F_Eout = H.*F_Ein; %multiply
%     Eout = ifftshift(ifft2(F_Eout)); %inv fft, center obs field
% end
% % 



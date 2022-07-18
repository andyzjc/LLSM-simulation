close all
clear all

% delta k = 1/wavelength

% Simulation of lattice light sheet pattern, see Chen 2014, S5
% Andy Zhang, 07/08/22
%% Initialize some parameters for back pupil plane
Nxz = 1024; % pixels
Ny = 64; % pixels

% Physical Parameter 
wavelength = 488 * 10^-6; % mm
n = 1.33; % water
NAdect = 1.1; 
NAmin = 0.5;
NAmax = 0.6;
NAideal = (NAmin + NAmax)/2;
f_obj = 10; % mm
D_obj = 12; % mm

% Real space image parameter
xz_scale = 4;
y_scale = 16;
deltax = wavelength / n / xz_scale; % 8 pixel / wavelength in medium 
deltaz = deltax;

% Fourier space image parameter
k0 = 2*pi/wavelength;
k0_med = k0 * n; 
deltakx = 2*pi / (Nxz * deltax);
deltakz = deltakx;
kmax = Nxz/2 * deltakx;

% Create mesh 
[ax, az] = meshgrid(  (-Nxz/2+1) : (Nxz/2) ) ; 
kx = deltakx * ax;
kz = deltakz * az;

x = deltax * ax;
z = deltaz * az; 

%% Generate Ideal lattice 
Illumi_ideal = zeros(size(kx));
theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [0, 90, 180, 270]; % square
% theta = [90, 270]; % standing wave
k_illum = round((kmax * NAideal)/deltakx); % pixels
k_NAmax = round((kmax * NAmax)/deltakx); % pixels
k_NAmin = round((kmax * NAmin)/deltakx); % pixels
k_diff = k_NAmax - k_NAmin; % pixels
kxposition = round(k_illum * cos(theta * pi/180));
kzposition = round(k_illum * sin(theta * pi/180));

% Ideal lattice
for j = 1:length(kxposition)
    Illumi_ideal(...
        (Nxz/2+1) + kzposition(j) - round(k_diff/2) : (Nxz/2+1) + kzposition(j) + round(k_diff/2),...
        (Nxz/2+1) + kxposition(j)) = 1;
end

%% Generate mask and 
A_mask = zeros(size(kx));
radius = k_illum; % pixels
thickness = k_diff;
A_mask = heaviside(thickness/2 * deltakx - abs( sqrt(kx.^2 + kz.^2 ) - radius * deltakx ));
I_mask = Illumi_ideal .* A_mask;
fig1 = figure(1);
    image1 = imagesc(kx(:,1), kz(1,:), I_mask);
    caxis([min(min(I_mask)) max(max(I_mask))]);
    colormap(jet)
    title("Rear Pupil Illumination")
    xlabel("kx (1/mm)")
    ylabel("ky (1/mm)")
    colorbar;
    axis image;


%% PSF(x,z)
PSF_exc = abs( fftshift(fft2(I_mask )) ).^2;
fig2 = figure(2);
    image2 = imagesc(x(:,1), z(1,:), PSF_exc);
    caxis([min(min(PSF_exc)) max(max(PSF_exc))]);
    colormap(jet)
    title("Excitation PSF")
    xlabel("x (mm)")
    ylabel("y (mm)")
    colorbar;
    axis image;

%% OTF
OTF_exc = abs(fftshift(fft2(PSF_exc)));
fig3 = figure(3);
    image3 = imagesc(kx(1,:),kz(:,1),OTF_exc);
    caxis([min(min(OTF_exc)) max(max(OTF_exc))]);
    colormap(jet)
    title("Excitation OTF")
    xlabel("kx (1/mm)")
    ylabel("ky (1/mm)")
    colorbar;
    axis image;

%% Dither 
PSF_exc_dither = meshgrid(sum(PSF_exc,2))';  

fig4 = figure(4);
    image5 = imagesc(x(1,:),z(:,1),PSF_exc_dither);
    caxis([min(min(PSF_exc_dither)) max(max(PSF_exc_dither))]);
    colormap(jet)
    title("Excitation PSF - Dither")
    xlabel("x (mm)")
    ylabel("y (mm)")
    colorbar;
    axis image;

%% Dither OTF 
OTF_dither = abs(fftshift(fft2(PSF_exc_dither)));

fig5 = figure(5);
    image7 = imagesc(kx(1,:),kz(:,1),OTF_dither);
    caxis([min(min(OTF_dither)) max(max(OTF_dither))]);
    colormap(jet)
    title("Excitation OTF - Dither")
    xlabel("kx (1/mm)")
    ylabel("ky (1/mm)")
    colorbar;
    axis image;

%% Propagation in y-direction
% Propagator
deltay = deltax * y_scale; % 1 pixel = 2 wavelength
ky = sqrt(k0_med^2 - kx.^2 - kz.^2);
y = (1:Ny) * deltay;
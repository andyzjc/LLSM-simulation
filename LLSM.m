close all
clear all

% delta k = 1/wavelength

% Simulation of lattice light sheet pattern, see Chen 2014, S5
% Andy Zhang, 07/08/22
%% Initialize some parameters for back pupil plane
Nxz = 512; % pixels
Ny = 64; % pixels

% Physical Parameter 
wavelength = 0.488; % um
n = 1.33; % water
NAdect = 1.1; 
NAmin = 0.3;
NAmax = 0.6;
NAideal = (NAmin + NAmax)/2;
f = 10 * 10^3; % mm
D = 12 * 10^3; % mm

% Real space image parameter
xz_scale = 4;
y_scale = 16;
deltax = wavelength / n / xz_scale; % 8 pixel / wavelength in medium 
deltaz = deltax;
deltay = deltax * y_scale; % 1 pixel = 2 wavelength

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
ky = sqrt(k0_med^2 - kx.^2 - kz.^2);
x = deltax * ax;
z = deltaz * az;
y = (1:Ny) * deltay; 

%% Generate Ideal lattice 
Illumi_ideal = zeros(size(kx));
theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [0, 90, 180, 270]; % square
% theta = [90, 270];
k_illum = 30; % pixels
kxposition = round(k_illum * cos(theta * pi/180));
kzposition = round(k_illum * sin(theta * pi/180));

for j = 1:length(kxposition)
    Illumi_ideal((Nxz/2+1) + kzposition(j), (Nxz/2+1) + kxposition(j)) = 1;
end
E_ideal = ifft2(fftshift(Illumi_ideal));
fig1 = figure(1);
image1 = imagesc(x(1,:),z(:,1),abs(E_ideal));
    title("Ideal Optical Lattice, E_{ideal}");
    xlabel("x (um)")
    ylabel("z (um)")
    colormap(jet);
    image1.Parent.YDir = 'normal';
    caxis([min(min(abs(E_ideal))), max(max(abs(E_ideal)))])
%     image1.Parent.XAxis.TickValues = [];
%     image1.Parent.YAxis.TickValues = [];
    colorbar

%%


% %% Calculate OTF
% OTF_exc = abs(fftshift(fft2(PSF_exc)));
% 
% fig4 = figure(4);
% image4 = imagesc(kx(1,:),kz(:,1),OTF_exc);
%     colormap(jet)
%     title("Excitation OTF")
%     xlabel("kx (1/um)")
%     ylabel("ky (1/um)")
% 
% %% Dither mode by swepting
% PSF_exc_dither = meshgrid(sum(PSF_exc,2))';  
% 
% fig5 = figure(5);
% image5 = imagesc(x(1,:),z(:,1),PSF_exc_dither);
%     colormap(jet)
%     title("Excitation PSF - Dither mode")
%     xlabel("x (um)")
%     ylabel("z (um)");
%     fig5.Children.YDir = 'normal';
% 
% %% Dither OTF 
% OTF_dither = abs(fftshift(fft2(PSF_exc_dither)));
% 
% fig7 = figure(7);
% image7 = imagesc(kx(1,:),kz(:,1),OTF_dither);
%     colormap(jet)
%     title("Excitation OTF, Dithered")
%     xlabel("kx (1/um)")
%     ylabel("ky (1/um)")



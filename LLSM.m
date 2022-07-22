close all
clear all

%% Initialize some parameters for back pupil plane
Nxz = 1025; % pixels
Ny = 257; % pixels

% Physical Parameter 
wavelength = 0.488; % um
n = 1.33; % water
NAdect = 1;
NAmin = 0.57;
NAmax = 0.65;
NAideal = (NAmin + NAmax)/2;
dither_period = 5; % in micron
dither_step = 201; % number of steps per dither period 

% when 2 plane wave traveling in opposite direction and interfere with each
% other, the periodicity of resulting standing wave is wavelength/4pi

% Real space image parameter
deltax = wavelength/4/pi; % oversampling
deltaz = deltax;

% Fourier space image parameter
k0_med = 2* pi/wavelength;
Kbound = 1/(2*deltax); % Nyquist
deltakx = 2*Kbound / Nxz;
deltakz = deltakx;
% k_ideal = NAideal / (wavelength); % in unit of 1/um
k_ideal = NAideal/ ( 2* pi) * Kbound;

% Create mesh 
[ax, az] = meshgrid(  -(Nxz-1)/2 : (Nxz-1)/2 ) ; 
kx = deltakx * ax;  %in unit 1/um
kz = deltakz * az;
x = deltax * ax;
z = deltaz * az; 
X = x(1,:);
Z = z(:,1);
KX = kx(1,:);
KZ = kz(:,1)';

%% Generate Ideal lattice 
Illumi_ideal = zeros(size(kx));
theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [30, 150, 210, 330]; % square
% theta = [90, 270]; % standing wave

kxposition = k_ideal * cosd(theta) / deltakx; % pixel
kzposition = k_ideal * sind(theta) / deltakz; % pixel

% kxposition_sign = sign(kxposition);
% kzposition_sign = sign(kzposition);
% 
% kxpositiontop = ceil(abs(kxposition)) .* kxposition_sign
% kxpositionbottom = floor(abs(kxposition)) .* kxposition_sign
% 
% kzpositiontop = ceil(abs(kzposition)) .* kzposition_sign
% kzpositionbottom = floor(abs(kzposition)) .* kzposition_sign
% 
% % 
% % kxposition_round = (Nxz-1)/2 + round(abs(kxposition)) .* kxposition_sign;
% % kzposition_round = (Nxz-1)/2 + round(abs(kzposition)) .* kzposition_sign;
% 
% kzposition = [78   157    78   -78  -157  -78]
% Ideal lattice illumination
for j = 1:length(kxposition)

    Illumi_ideal( ...
        (Nxz-1)/2 + round(kzposition(j)) ,...
        (Nxz-1)/2 + round(kxposition(j)) ) = 1;
end

% Ideal lattice
E_ideal = myFFT(Illumi_ideal,0);

fig1 = figure(1);
    subplot(1,2,1)
image100 = imagesc(KX, KZ, Illumi_ideal);
    colormap("gray")
    title("Rear Pupil Illumination")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

    subplot(1,2,2)
image1 = imagesc(X, Z ,abs(E_ideal));
    caxis([min(min(abs(E_ideal))) max(max(abs(E_ideal)))]);
    colormap("jet")
    title("Ideal Lattice")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

%% Gaussian Bounding and get rear pupil illum back
a = 4*wavelength;
gauss_bound = exp(-2 * z.^2 / a^2);
E_bound = gauss_bound .* E_ideal;
fig2 = figure(2);
    subplot(1,2,1)
    image2 = imagesc(X, Z,abs(E_bound));
    colormap("jet")
    title("Bounded Ideal Lattice")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

Illum_bound = abs(myFFT(E_bound,1)).^2;
Illum_bound = Illum_bound/max(max(Illum_bound));
    subplot(1,2,2)
    image3 = imagesc(KX, KZ,Illum_bound);
    colormap(jet)
%     caxis([min(min(Illum_bound)) max(max(Illum_bound))]);
    title("Bounded Illumination")
    xlabel("kx (Constant/lambda)")
    ylabel("kz (Constant/lambda)")
    colorbar;
    axis image;

%% Generate mask and filter
A_mask = zeros(size(kx));
k_NAmax = NAmax/wavelength; % 1/um
k_NAmin = NAmin/wavelength; 

A_mask = ((k_NAmax > sqrt(kx.^2 + kz.^2)) .* (k_NAmin < sqrt(kx.^2 + kz.^2)));
fig4 = figure(4);
    subplot(1,2,1)
    image4 = imagesc(KX, KZ, imfuse(A_mask,Illum_bound));
    title("Mask")
    xlabel("kx (Constant/lambda)")
    ylabel("kz (Constant/lambda)")
    colorbar;
    axis image;
Pupil_fun = Illum_bound .* A_mask;
    subplot(1,2,2)
    image1 = imagesc(KX, KZ, Pupil_fun);
    caxis([min(min(Pupil_fun)) max(max(Pupil_fun))]);
    colormap(jet)
    title("Pupil Function")
    xlabel("kx (Constant/lambda)")
    ylabel("kz (Constant/lambda)")
    colorbar;
    axis image;


%% PSF(x,z)
% PSF_exc = abs( fftshift(fft2(fftshift(Pupil_fun) )) ).^2;
PSF_exc = abs( fftshift(fft2(fftshift(Pupil_fun))) ).^2;
PSF_exc = PSF_exc/max(max(PSF_exc));
fig5 = figure(6);
    image2 = imagesc(X, Z,PSF_exc);
    caxis([min(min(PSF_exc)) max(max(PSF_exc))]);
    colormap(jet)
    title("Excitation PSF")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;
    axis image;

%% OTF
OTF_exc = abs(fftshift(fft2(fftshift(PSF_exc))));
OTF_exc = OTF_exc/max(max(OTF_exc));
fig6 = figure(6);
    image3 = imagesc(KX, KZ,OTF_exc);
    caxis([min(min(OTF_exc)) max(max(OTF_exc))]);
    colormap(jet)
    title("Excitation OTF")
    xlabel("kx (Constant/lambda)")
    ylabel("kz (Constant/lambda)")
    colorbar;
    axis image;

%% Dither 
dither_period = 10/deltax; % micron
PSF_exc_dither = PSF_exc;
for jj = 1:dither_step
    PSF_exc_dither = PSF_exc_dither + circshift(PSF_exc,round(jj*dither_period/dither_step),2)/2;
end
% PSF_exc_dither = meshgrid(sum(PSF_exc(:,1:dither_period/deltax),2))';  
PSF_exc_dither = PSF_exc_dither / max(max(PSF_exc_dither));

fig7 = figure(7);

    subplot(1,2,2)
h1 = plot(Z,PSF_exc_dither(:,(Nxz-1)/2));
    h1.Parent.YAxis.TickValues = [0, 0.25, 0.5, 0.75, 1];
    grid on

subplot(1,2,1)
image4 = imagesc(X, Z ,PSF_exc_dither);
    colormap(jet)
    title("Excitation PSF - Dither")
    xlabel("x/wavelength ")
    ylabel("z/wavelength ")
    colorbar;length(dither_step)
    axis image;

%% Dither OTF 
OTF_dither = abs(fftshift(fft2(PSF_exc_dither)));

fig8 = figure(8);
    image5 = imagesc(kx(1,:)/Kbound, kz(:,1)/Kbound, OTF_dither);
    caxis([min(min(OTF_dither)) max(max(OTF_dither))]);
    colormap(jet)
    title("Excitation OTF - Dither")
    xlabel("kx (Constant/lambda)")
    ylabel("kz (Constant/lambda)")
    colorbar;
    axis image;

%% Propagation in y-direction
y_scale = 100; 
deltay = wavelength/4/pi * y_scale;
ky = sqrt(k0_med^2 - kx.^2 - kz.^2);
y = (0:Ny) * deltay; % max(y) = FOVy 

E_prop = zeros(Nxz,Nxz,Ny);
I_prop = E_prop;
Ein = Pupil_fun/(max(max(Pupil_fun)));

%% Propagator
% for i = 1:length(y)
for i = 1:length(y)
    propagator = exp(1i * ky * y(i));
    Eout = fftshift(fft2((Ein .* propagator)));
    imagesc(x(1,:)/wavelength,z(:,1)/wavelength,abs(Eout.^2))
    drawnow
    colormap(jet)
    pause(0.5)
    caxis([0 4000])
    disp(i)
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



clear all
close all

%% Physical parameters
 theta = [30, 90, 150, 210, 270, 330]; % hexogonal
% theta = [30, 150, 210, 330]; % square
% theta = [90, 270]; % standing wave

weighting =  [1, 1, 1, 1, 1, 1]; % beam weighting

% Physical Parameter 
N = 513; % pixels
n = 1.33;
lambda_exc = 0.488; % um 
lambda_det = 0.488;
wavelength_exc =  lambda_exc / n;
wavelength_det = lambda_det / n;
NAmin = 0.57;
NAmax = 0.65;
NAdet = 1.1;
NAideal = (NAmin + NAmax)/2;
dither_period = 3; % um
dither_step = 201; % number of s teps per dither period 
gauss_bound_width = 1; % Gaussian Bounding, um
xz_scale = 4;
y_scale = 2;

k_wave = 1/wavelength_exc;
k_ideal = k_wave * NAideal / n;
k_det = k_wave * NAdet / n;
k_bound = xz_scale * k_wave;
k_NAmax = NAmax /n * k_wave; % k
k_NAmin = NAmin /n * k_wave; 
deltak = 2 * k_bound / N;
deltax = 1/(2 * k_bound);

% excitation
[ax, az] = meshgrid(  -(N-1)/2 : (N-1)/2 ) ; 
kx_exc = deltak * ax;  %in unit wavelength
kz_exc = kx_exc';
ky_exc = sqrt(k_wave^2 - kx_exc.^2 - kz_exc.^2);
ky_exc(kx_exc.^2 + kz_exc.^2 > k_wave.^2 ) = 0;
x_exc = deltax * ax; 
z_exc = x_exc'; 
y_exc = (-(N+1)/2+1 : (N+1)/2-1) * deltax * y_scale; 

% detection
kx_det = kx_exc;
ky_det = kz_exc;
kz_det = ky_exc;
x_det = x_exc;
y_det = z_exc;
z_det = y_exc;

% for displaying
KX_exc = kx_exc(1,:) / k_wave;
KZ_exc = KX_exc';
X_exc = x_exc(1,:)  / wavelength_exc; % value * wavelength = physical value (um)
Z_exc = X_exc'; 
Y_exc = y_exc / wavelength_exc; 

KX_det = KX_exc;
KY_det = KZ_exc;
X_det = X_exc;
Y_det = Z_exc;
Z_det = Y_exc;

%% Simulation 
% Generate Ideal lattice back pupil 
Illumi_ideal = zeros(size(ax));
kxposition = k_ideal * cosd(theta) /deltak; % pixel
kzposition = k_ideal * sind(theta) /deltak; % pixel

for j = 1:length(kxposition)

    Illumi_ideal( ...
        (N+1)/2 + round(kzposition(j)) ,...
        (N+1)/2 + round(kxposition(j)) ) = 1 * weighting(j);
end
E_ideal = ifft2(Illumi_ideal); 

% bounded lattice 
gauss_bound = exp(-2 * z_exc.^2 / (gauss_bound_width)^2);
E_bound = gauss_bound .* E_ideal;

% bounded back pupil
Illum_bound = abs(fft2(fftshift(E_bound))).^2;
Illum_bound = Illum_bound/max(max(Illum_bound));

% Generate mask 
A_mask = ((k_NAmax > sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_NAmin < sqrt(kx_exc.^2 + kz_exc.^2)));

% Pupil functions
Pupil_fun_exc = Illum_bound .* A_mask;
Pupil_fun_det = k_det > sqrt(kx_det.^2 + ky_det.^2);

% % Propagation
tic
PSF_exc_3d = zeros(N,N, N);
OTF_exc_3d = zeros(N,N, N);
OTF_exc_3d_phase = zeros(N,N, N);
PSF_exc_3d_dither = zeros(N,N, N);
OTF_exc_3d_dither = zeros(N,N, N);
OTF_exc_3d_dither_phase = zeros(N,N, N);
PSF_det_3d = zeros(N,N,N);

% propagation
for i = 1:length(y_exc)
    propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
    PSF_exc_3d(:,:,i) = abs( fftshift( ifft2(Pupil_fun_exc .* propagator_exc) ) ).^2;
    OTF_exc_3d(:,:,i) = abs(fftshift(fft2(PSF_exc_3d(:,:,i))));
    OTF_exc_3d_phase(:,:,i) = angle(OTF_exc_3d(:,:,i));
end  

% detection propagation
for ii = 1:length(z_det)
    propagator_det = exp(2*pi * 1i * kz_det * z_det(ii));
    PSF_det_3d(:,:,ii) = abs( fftshift( ifft2(Pupil_fun_det .* propagator_det) ) ).^2;
end

% dithering along x exc
for j = 1:dither_step
    PSF_exc_3d_dither = PSF_exc_3d_dither + ...
        circshift(PSF_exc_3d,round(j * dither_period / deltax / dither_step),2);
end

for k = 1:length(y_exc)
    OTF_exc_3d_dither(:,:,k) = abs(fftshift(fft2(PSF_exc_3d_dither(:,:,k))));
%     OTF_exc_3d_dither_phase(:,:,k) = angle( OTF_exc_3d_dither(:,:,k) );
end

% Overall 
Overall_PSF_axial = squeeze(PSF_exc_3d(:,:,(N+1)/2)) .* squeeze(PSF_det_3d(:,(N+1)/2,:))' ; 
Overall_PSF_lateral = squeeze(PSF_exc_3d(:,(N+1)/2,:)) .* squeeze(PSF_det_3d(:,:,(N+1)/2));
Overall_OTF_axial = abs(fftshift(fft2(Overall_PSF_axial)));
Overall_OTF_lateral =  abs(fftshift(fft2(Overall_PSF_lateral)));

% Normalize
PSF_exc_3d = PSF_exc_3d/max(max(max(PSF_exc_3d)));
OTF_exc_3d = OTF_exc_3d/max(max(max(OTF_exc_3d)));
PSF_exc_3d_dither = PSF_exc_3d_dither/max(max(max(PSF_exc_3d_dither)));
OTF_exc_3d_dither = OTF_exc_3d_dither/max(max(max(OTF_exc_3d_dither)));
PSF_det_3d = PSF_det_3d/max(max(max(PSF_det_3d)));
Overall_PSF_axial = Overall_PSF_axial/max(max(Overall_PSF_axial));
Overall_PSF_lateral = Overall_PSF_lateral/max(max(Overall_PSF_lateral));
Overall_OTF_axial = Overall_OTF_axial/max(max(Overall_OTF_axial));
Overall_OTF_lateral = Overall_OTF_lateral/max(max(Overall_OTF_lateral));
toc

%% Figure 1 - Rear Pupil 
    fig1 = figure(1);
    fig1.Name = "XZ-excitation, Y = 0";
    colormap(hot(256))

     subplot(3,4,1);
image11 = imagesc(KX_exc,KZ_exc, real(Illumi_ideal) );
    title("Ideal Lattice, " +...
          "N_{xz} = " + num2str(N) + ...
          ", K_{bound} = " + num2str(k_bound) + " (1/um)")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image11.Parent.XLim = [-1,1];
    image11.Parent.YLim = [-1,1];
    colorbar;    

    subplot(3,4,2);
image12 = imagesc(X_exc,Z_exc, abs( E_ideal));
    title("Ideal Lattice," + ...
          "\lambda_{exc}/n = " + num2str(lambda_exc, '%.3f') + "um / " + ...
           num2str(n))
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;

    subplot(3,4,3);
image13 = imagesc(X_exc, Z_exc, abs(E_bound));
    title("Bounded Ideal Lattice, " + ...
          "Bound width = " + num2str(gauss_bound_width) + "um")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;
    
    subplot(3,4,4)
image14 = imagesc( KX_exc, KZ_exc,...
                  Illum_bound );
    title("Rear pupil illumination after bounding")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image14.Parent.XLim = [-1,1];
    image14.Parent.YLim = [-1,1];
    colorbar;

Illum_mask = imfuse(Illum_bound,A_mask,"falsecolor","ColorChannels","green-magenta");
    subplot(3,4,5)
image15 = imagesc( KX_exc, KZ_exc,...
                  Illum_mask);
    title("Masking, " +...
          "NA_{max} = " + num2str(NAmax) +...
          ", NA_{min} = " + num2str(NAmin) )
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image15.Parent.XLim = [-1,1];
    image15.Parent.YLim = [-1,1];

    subplot(3,4,6)
image16 = imagesc( KX_exc, KZ_exc,...
                  Pupil_fun_exc );
    title("Bounded Rear Pupil")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image16.Parent.XLim = [-1,1];
    image16.Parent.YLim = [-1,1];
    colorbar;

    subplot(3,4,7)
image17 = imagesc(X_exc, Z_exc, PSF_exc_3d(:,:,(N+1)/2));
    title("XZ-Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(3,4,8)
image18 = imagesc( KX_exc,...
                  KZ_exc,...
                  OTF_exc_3d(:,:,(N+1)/2) ) ;
    title("XZ-Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image
    image18.Parent.XLim = [-1,1];
    image18.Parent.YLim = [-1,1];


    subplot(3,4,9)
    hold on;
image19 = imagesc(X_exc, Z_exc ,PSF_exc_3d_dither(:,:,(N+1)/2));
    title("Dithered XZ-Excitation PSF, " + ...
          "T_d = " + num2str(dither_period) + "um, " +...
          "Dither Step = " + num2str(dither_step))
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(3,4,10)
    zPSF = squeeze(PSF_exc_3d_dither(:,(N+1)/2,(N+1)/2))/max(squeeze(PSF_exc_3d_dither(:,(N+1)/2,(N+1)/2)));
image110 = plot( Z_exc, zPSF);
    title("Dithered Z-Excitation PSF")
    ylabel("Normalized a.u. ")
    xlabel("z/\lambda")
    image110.LineWidth = 2;
    image110.Color = 'r';
    image110.Parent.YAxis.TickValues = linspace(0,1,11);
    image110.Parent.XLim = [-15,15];
    grid on

    subplot(3,4,11)
image111 = imagesc( KX_exc,...
                   KZ_exc,...
                   OTF_exc_3d_dither(:,:,(N+1)/2)  );
    title("Dithered XZ-Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image
    image111.Parent.XLim = [-1,1];
    image111.Parent.YLim = [-1,1];

    subplot(3,4,12)
image112 = plot( KZ_exc, squeeze(OTF_exc_3d_dither(:,(N+1)/2,(N+1)/2)));
    title("Dithered XZ-Excitation OTF, " + ...
        "K_x = " + num2str(0) )
    ylabel("Normalized a.u. ")
    xlabel("kz * \lambda")
    image112.LineWidth = 2;
    image112.Color = 'r';
    image112.Parent.XLim = [-2,2];
    image112.Parent.YAxis.TickValues = linspace(0,1,11);
    grid on

%% Figure 2 - Excitation
    fig2 = figure(2);
    fig2.Name = "Focal PSF/OTF";
    colormap(hot(256))

     subplot(2,3,1);
image21 = imagesc(X_exc,Z_exc,PSF_exc_3d(:,:,(N+1)/2) );
    title("XZ-Excitation PSF, "  + "Y = 0")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(2,3,2);
image22 = imagesc(X_exc,Z_exc,PSF_exc_3d_dither(:,:,(N+1)/2));
     title("Dithered XZ-Excitation PSF, " + ...
          "T_d = " + num2str(dither_period) + "um, " +...
          "Step = " + num2str(dither_step) + ...
          ", Y = 0")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;
    
    subplot(2,3,3);
image23 = plot( KZ_exc, OTF_exc_3d(:,(N+1)/2,(N+1)/2));
    title("Z-Excitation-OTF, " + "K_X = 0, " + "K_Y = 0")
    ylabel("Normalized a.u. ")
    xlabel("kz * \lambda")
    image23.Color = 'r';
    image23.LineWidth = 2;
    image23.Parent.XLim = [-2,2];
    colorbar;
    grid on

    subplot(2,3,4);
image24 = imagesc(Y_exc, Z_exc, squeeze(PSF_exc_3d(:,(N+1)/2,:)));
    title("YZ-Excitation-PSF, " + "X = 0" )
    axis image
    xlabel("y/\lambda")
    ylabel("z/\lambda")
    colorbar

    subplot(2,3,5);
image25 = imagesc(Y_exc,Z_exc, squeeze(PSF_exc_3d_dither(:,(N+1)/2,:)) );
    title("Dithered YZ-Excitation-PSF, " + "X = 0" )
    xlabel("y/\lambda")
    ylabel("z/\lambda")
    axis image
    colorbar;

    subplot(2,3,6);
    yPSF_exc = squeeze(PSF_exc_3d((N+1)/2,(N+1)/2,:));
    % Calculate yFWHM
    index = find(yPSF_exc >= 0.5);
    yFWHM = Y_exc(index(end)) - Y_exc(index(1)); % pixels
image26 = plot(Y_exc, yPSF_exc );
    title("Y-Excitation-PSF, " + "X = 0, Z = 0, " + ...
          "yFWHM = " + num2str(yFWHM) + "\lambda")
    xlabel("y/\lambda")
    ylabel("Normalized a.u. ")
    image26.Color = 'r';
    image26.LineWidth = 2;
    colorbar;
    grid on
    
%% Figure 3 - Overall PSF/OTF
    fig3 = figure(3);  
    fig3.Name = "Overall XZ-Axial PSF/OTF, focal plane";
    colormap(hot(256))
        
    subplot(2,4,1)
image31 = imagesc(X_exc,Z_exc,PSF_exc_3d_dither(:,:,(N+1)/2) );
    title("Dithered XZ-Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(2,4,2)
 image32 = imagesc(X_det,Z_det,squeeze(PSF_det_3d(:,(N+1)/2,:))');
    title("XZ-Detection PSF ")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;  
    image32.Parent.XLim = [-5,5];
    image32.Parent.YLim = [-5,5];

    subplot(2,4,3)
 image33 = imagesc(X_exc,Z_exc,Overall_PSF_axial);
    title("Overall PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;  
    image33.Parent.XLim = [-5,5];
    image33.Parent.YLim = [-5,5];

    subplot(2,4,5)
image35 = imagesc(KX_exc,...
                  KZ_exc,...
                  OTF_exc_3d(:,:,(N+1)/2) );
    title("XZ-Excitation OTF ")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;
    image35.Parent.XLim = [-2,2];
    image35.Parent.YLim = [-2,2];

    subplot(2,4,6)
 image36 = imagesc(KX_det,...
                  KZ_exc,...
                  abs(fftshift(fft2(squeeze(PSF_det_3d(:,(N+1)/2,:))'))));
    title("XZ-Detection OTF ")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;  
    image36.Parent.XLim = [-2,2];
    image36.Parent.YLim = [-2,2];

    subplot(2,4,7)
 image37 = imagesc(X_exc,Z_exc,Overall_OTF_axial);
    title("Overall OTF")
    xlabel("kx * \la5mbda")
    ylabel("kz * \lambda")
    colorbar;
    axis image;  
    image37.Parent.XLim = [-20,20];
    image37.Parent.YLim = [-20,20];

    h1 = subplot(2,4,[4,8]);
    hold on
line_exc = plot(zPSF,Z_exc);
    line_exc.Color = 'g';
    line_exc.LineWidth = 2;
line_det = plot(squeeze(PSF_det_3d((N+1)/2,(N+1)/2,:)),Z_exc);
    line_det.Color = 'b';
    line_det.LineWidth = 2;
line_overall = plot(squeeze(Overall_PSF_axial(:,(N+1)/2)),Z_exc);
    line_overall.Color = 'r';
    line_overall.LineWidth = 2;
% line_lateral = plot(squeeze(PSF_det_3d(:,(N+1)/2,(N+1)/2)),Z_exc);
%     line_lateral.Color = 'k';
%     line_lateral.LineWidth = 2;
    title("Overall Axial-PSF")
    ylabel("z/\lambda")
    xlabel("Normalized a.u. ")
    lgd = legend("Excitation", "Detection","Overall");  
        lgd.FontWeight = 'bold';
        lgd.FontSize = 7;
        lgd.LineWidth = 1;
    h1.YLim = [-6,6];
    h1.YTick = linspace(-6,6,13);
    grid on
    hold off

    

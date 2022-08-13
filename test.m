clear all
% Physical Parameter 
N = 257; % pixels
n = 1.33;
tic
lambda_det = 0.515;
wavelength_det = lambda_det / n;
NAdet = 1.0;
k_xz_scale = 8;
k_wave = 1/wavelength_det  ;
k_det = k_wave * NAdet / n ;
k_bound = k_xz_scale * k_wave;
deltak = 2 * k_bound / N;
deltax = 1/(2 * k_bound);

% detection
 [ax, ay] = meshgrid(  -(N-1)/2 : (N-1)/2 ) ; 
kx_det = deltak * ax;  %in unit wavelength
ky_det = kx_det';

kz_det = sqrt(k_wave^2 - kx_det.^2 - ky_det.^2) .* (k_det.^2 >= kx_det.^2 + ky_det.^2);
kz_det(kx_det.^2 + ky_det.^2 >= k_wave.^2 ) = 0;

x_det = deltax * ax; 
y_det = x_det'; 
z_det =  (  -(N-1)/2 : (N-1)/2 ) * deltax;
KZ_det = (  -(N-1)/2 : (N-1)/2 ) * 1/(2*max(z_det)) / k_wave;

% for displaying
KX_det = kx_det(1,:) / k_wave;
KY_det = KX_det'; 
X_det = x_det(1,:)  / (wavelength_det); 
Y_det = X_det'; 
Z_det = z_det / (wavelength_det); 

% Pupil_fun_det = double(k_det.^2 >= kx_det.^2 + ky_det.^2);
Pupil_fun_det = k_wave./kz_det;
Pupil_fun_det(Pupil_fun_det == Inf) = 0;

%%
% detection propagation
PSF_det_3d = zeros(N,N,N);
for ii = 1:length(z_det) 
    propagator_det = exp(2*pi * 1i * kz_det * z_det(ii)); 
    PSF_det_3d(:,:,ii) = abs( fftshift( ifft2(Pupil_fun_det .* propagator_det) ) ).^2;  
end 

% halfpix = (N-1)/2;
% halfrange = ceil(k_det/deltak);
% NumPointsWithinBoundingSquare = (2.*halfrange + 1).^2;
% %define matrices to store the info about the illuminated points in the pupil:
% NAx = zeros(1,NumPointsWithinBoundingSquare);
% NAy = NAx;
% PupilEfield = NAx;
% numpoints = 1;
% EAmpThreshhold = 0.001;
% for i = (halfpix-halfrange):(halfpix+halfrange)
%     for j = (halfpix-halfrange):(halfpix+halfrange)
%         if Pupil_fun_det(i,j) > EAmpThreshhold
%             NAx(numpoints) = i;
%             NAy(numpoints) = j;
%             PupilEfield(numpoints) = Pupil_fun_det(i,j);
%             numpoints = numpoints + 1;
%             kxindex(numpoints) = kx_det(i,j);
%             kyindex(numpoints) = ky_det(i,j);
%             kzindex(numpoints) = kz_det(i,j);
%         end
%     end
% end
% numpoints = numpoints-1;
% NAx = NAx(1:numpoints);
% NAy = NAy(1:numpoints);
% 
% for ii = 1:length(z_det) 
%     Efield = zeros(N,N);
%     
%     for jj = 1:numpoints
%         Efield = Efield + PupilEfield(jj) .* exp(2.*pi.*1i .* (kxindex(jj).*x_det + kyindex(jj).* y_det + kzindex(jj).*z_det(ii)));
%     end
%     PSF_det_3d_integration(:,:,ii) = abs(Efield).^2;
% end 

% NAx = KX_det(NAx) * n / k_wave;
% NAy = KY_det(NAy) * n / k_wave;
% sinth = sqrt(NAx.^2 + (NAy').^2);
% costh = sqrt(1 - sinth.^2);  %angle of k vector to y axis
% sinphi = NAy./sinth;
% cosphi = NAx./sinth;
% kx = cosphi.*sinth;
% kz = sinphi.*sinth;
% ky = costh;
% 
% PupilEfield = PupilEfield(1:numpoints);
% ax = linspace(-50,50,N)
% az = ax';
% az = (  -(N-1)/2 : (N-1)/2 );
% for ii = 1:length(z_det) 
%     Efield = zeros(N,N);
%     
%     for jj = 1:numpoints
%         Efield = Efield + PupilEfield(jj) .* exp(2.*pi.*1i .* (kx(jj).*ax + kz(jj).* ay + ky(jj).*az(ii)));
%     end
%     PSF_det_3d_integration(:,:,ii) = abs(Efield).^2;
% end 

% PSF_det_3d = PSF_det_3d./max(max(max(PSF_det_3d)));
%%
fig1 = figure(1);
    colormap(hot)
    subplot(2,3,1)
xyPSF1 = squeeze(PSF_det_3d(:,:,(N+1)/2));
xyPSF1 = xyPSF1/max(max(xyPSF1));
imagesc(X_det,Y_det,xyPSF1);
    title("Lateral PSF")
    xlabel("x / \lambda/n")
    ylabel("y / \lambda/n")
    axis image
    colorbar
    
    subplot(2,3,2)
xzPSF1 = squeeze(PSF_det_3d(:,(N+1)/2,:))';
xzPSF1 = xzPSF1/max(max(xzPSF1));
imagesc(X_det,Z_det,xzPSF1);
    title("Axial PSF")
    xlabel("x / \lambda/n")
    ylabel("z / \lambda/n")
    axis image
    colorbar
% 
    subplot(2,3,3)
xyOTF1 = abs(fftshift(fft2(xyPSF1)));
xyOTF1 = xyOTF1/max(max(xyOTF1));
imagesc(KX_det,KY_det,xyOTF1);
    title("Lateral OTF")
    xlabel("kx / (n/\lambda)")
    ylabel("ky / (n/\lambda)")  
    axis image
    colorbar

    subplot(2,3,4)
OTF_3d1 = abs(fftshift(fftn(PSF_det_3d)));
OTF_3d1(257,257,257) = 0;
OTF_3d1 = OTF_3d1/max(max(max(OTF_3d1)));
OTF = squeeze(OTF_3d1(:,(N+1)/2,:)); 
h0 = imagesc(KX_det,KZ_det,OTF'); 
    title("Axial OTF by 3d fftn")
    axis image
    h0.Parent.XLim = [-1,1];
    h0.Parent.YLim = [-1,1];
    colorbar

    subplot(2,3,5)
xzPSF_sum = squeeze(sum(PSF_det_3d(:,1:end,:)));
xzOTF1 = abs(fftshift(fft2(xzPSF_sum)))';
h1 = imagesc(KX_det,KZ_det,xzOTF1); 
    title("Axial OTF by 2d fft over sum of all x")
    xlabel("kx / (n/\lambda)") 
    ylabel("kz / (n/\lambda)")
    axis image
    h1.Parent.XLim = [-1,1];
    h1.Parent.YLim = [-1,1];
    colorbar

        subplot(2,3,6)
xzOTF2 = abs(fftshift(fft2(xzPSF1)))';
h2 = imagesc(KX_det,KZ_det,xzOTF2'); 
    title("Axial OTF by 2d fft over x = 0" )
    xlabel("kx / (n/\lambda)") 
    ylabel("kz / (n/\lambda)")
    axis image
    h2.Parent.XLim = [-1,1];
    h2.Parent.YLim = [-1,1];
    colorbar
  
    %%
%     fig1 = figure(1);
%     colormap(hot(256))
%     subplot(2,5,1)
%     xyPSF1 = squeeze(PSF_det_3d(:,:,(N+1)/2));
%     xyPSF1 = xyPSF1/max(max(xyPSF1));
% imagesc(X_det,Y_det,xyPSF1);
%     title("Lateral PSF")
%     xlabel("x / \lambda/n")
%     ylabel("y / \lambda/n")
%     axis image
%     colorbar
% 
%     subplot(2,5,2)
%     xzPSF1 = squeeze(PSF_det_3d(:,(N+1)/2,:))';
%     xzPSF1 = xzPSF1/max(max(xzPSF1));
% imagesc(X_det,Z_det,xzPSF1);
%     title("Axial PSF")
%     xlabel("x / \lambda/n")
%     ylabel("z / \lambda/n")
%     axis image
%     colorbar
% 
%     subplot(2,5,3)
% xyOTF1 = abs(fftshift(fft2(xyPSF1)));
% xyOTF1 = xyOTF1/max(max(xyOTF1));
% imagesc(KX_det,KY_det,xyOTF1);
%     title("Lateral OTF")
%     xlabel("kx / (n/\lambda)")
%     ylabel("ky / (n/\lambda)")  
%     axis image
%     colorbar
% 
%     subplot(2,5,4)
% OTF_3d1 = abs(fftshift(fftn(PSF_det_3d)));
% OTF_3d1 = OTF_3d1/max(max(max(OTF_3d1)));
% OTF = squeeze(OTF_3d1(:,(N+1)/2,:)); 
% imagesc(OTF); 
%     title("Axial OTF by 3d fft")
%     xlabel("kx / (n/\lambda)") 
%     ylabel("kz / (n/\lambda)")    
%     colorbar
%     axis image
% 
%     subplot(2,5,5)
% xzOTF1 = abs(fftshift(fft2(xzPSF1)));
% h1 = imagesc(KX_det,KZ_det,xzOTF1); 
%     title("Axial OTF by 2d fft")
%     xlabel("kx / (n/\lambda)") 
%     ylabel("kz / (n/\lambda)")
%     axis image
%     h1.Parent.XLim = [-1,1];
%     h1.Parent.YLim = [-1,1];
%     colorbar
% %
% subplot(2,5,6)
%     xyPSF2 = squeeze(PSF_det_3d_integration(:,:,(N+1)/2));
%     xyPSF2 = xyPSF2/max(max(xyPSF2));
% imagesc(X_det,Y_det,xyPSF2);
%     title("Lateral PSF")
%     xlabel("x / \lambda/n")
%     ylabel("y / \lambda/n")
%     axis image
%     colorbar
% 
%     subplot(2,5,7)
%     xzPSF2 = squeeze(PSF_det_3d_integration(:,(N+1)/2,:))';
%     xzPSF2 = xzPSF2/max(max(xzPSF2));
% imagesc(X_det,Z_det,xzPSF2);
%     title("Axial PSF")
%     xlabel("x / \lambda/n")
%     ylabel("z / \lambda/n")
%     axis image
%     colorbar
% 
%     subplot(2,5,8)
% xyOTF2 = abs(fftshift(fft2(xyPSF2)));
% xyOTF2 = xyOTF2/max(max(xyOTF2));
% imagesc(KX_det,KY_det,xyOTF2);
%     title("Lateral OTF")
%     xlabel("kx / (n/\lambda)")
%     ylabel("ky / (n/\lambda)")  
%     axis image
%     colorbar
% 
%     subplot(2,5,9)
% OTF_3d2 = abs(fftshift(fftn(PSF_det_3d_integration)));
% OTF_3d2 = OTF_3d2/max(max(max(OTF_3d2)));
% OTF = squeeze(OTF_3d2(:,(N+1)/2,:)); 
% imagesc(OTF); 
%     title("Axial OTF by 3d fft")
%     xlabel("kx / (n/\lambda)") 
%     ylabel("kz / (n/\lambda)")    
%     colorbar
%     axis image
% 
%     subplot(2,5,10)
% xzOTF2 = abs(fftshift(fft2(xzPSF2)));
% h1 = imagesc(KX_det,KZ_det,xzOTF2); 
%     title("Axial OTF by 2d fft")
%     xlabel("kx / (n/\lambda)") 
%     ylabel("kz / (n/\lambda)")
%     axis image
%     h1.Parent.XLim = [-1,1];
%     h1.Parent.YLim = [-1,1];
%     colorbar

%     subplot(2,5,5)
% xzOTF1 = abs(fftshift(fft2(xzPSF1,N*2-1,N*2-1)));
% h1 = imagesc(xzOTF1); 
%     title("Axial OTF")
%     xlabel("kx / (n/\lambda)") 
%     ylabel("kz / (n/\lambda)")
%     axis image
%     h1.Parent.XLim = [-1,1];
%     h1.Parent.YLim = [-1,1];
%     colorbar
%     
% subplot(2,5,6)
%     xyPSF2 = squeeze(PSF_det_3d_integration(:,:,(N+1)/2));
%     xyPSF2 = xyPSF2/max(max(xyPSF2));
% imagesc(xyPSF2);
%     title("Lateral PSF")
%     xlabel("x / \lambda/n")
%     ylabel("y / \lambda/n")
%     axis image
%     colorbar
% 
%     subplot(2,5,7)
%     xzPSF2 = squeeze(PSF_det_3d_integration(:,(N+1)/2,:))';
%     xzPSF2 = xzPSF2/max(max(xzPSF2));
% imagesc(xzPSF2);
%     title("Axial PSF")
%     xlabel("x / \lambda/n")
%     ylabel("z / \lambda/n")
%     axis image
%     colorbar
% % 
%     subplot(2,5,8)
% xyOTF2 = abs(fftshift(fft2(xyPSF2)));
% xyOTF2 = xyOTF2/max(max(xyOTF2));
% imagesc(xyOTF2);
%     title("Lateral OTF")
%     xlabel("kx / (n/\lambda)")
%     ylabel("ky / (n/\lambda)")  
%     axis image
%     colorbar
% 
%     subplot(2,5,9)
% OTF_3d2 = abs(fftshift(fftn(PSF_det_3d_integration)));
% OTF_3d2 = OTF_3d2/max(max(max(OTF_3d2)));
% OTF = squeeze(OTF_3d2(:,(N+1)/2,:)); 
% imagesc(OTF); 
%     title("Axial OTF")
%     xlabel("kx / (n/\lambda)") 
%     ylabel("kz / (n/\lambda)")    
%     colorbar
%     axis image
% 
%     subplot(2,5,10)
% xzOTF2 = abs(fftshift(fft2(xzPSF2,N*2-1,N*2-1)));
% h1 = imagesc(xzOTF2); 
%     title("Axial OTF")
%     xlabel("kx / (n/\lambda)") 
%     ylabel("kz / (n/\lambda)")
%     axis image
%     h1.Parent.XLim = [-1,1];
%     h1.Parent.YLim = [-1,1];
%     colorbar
%     
% PSF_det_3d_integration = PSF_det_3d_integration/max(max(max(PSF_det_3d_integration)));
PSF_det_3d = PSF_det_3d/max(max(max(PSF_det_3d)));
toc





function PrettyPlots
    global Data
    global Lattice
    global Figures

    KX_exc = Data.KX_exc;
    KZ_exc = Data.KZ_exc;
    KY_exc = Data.KY_exc;
    KX_det = Data.KX_det;
    KY_det = Data.KY_det;
    KZ_det = Data.KZ_det;
    X_exc = Data.X_exc;
    Z_exc = Data.Z_exc;
    Y_exc = Data.Y_exc;
    X_det = Data.X_det;
    Z_det = Data.Z_det;
    N = Data.N;
    n = Data.n;

    if Lattice.islattice == 1
        PlotLattice;
    else
        PlotOther;
    end

%% Figure 2 - Excitation
    Figures.fig2 = figure(2);
    Figures.fig2.Name = "yz-Excitation";
    Figures.fig2.WindowState = 'maximized';
    colormap(hot(256)) 

     subplot(2,3,1);
xzPSF_exc = Data.PSF_exc_3d(:,:,(Data.N+1)/2);
image21 = imagesc(X_exc,Z_exc,xzPSF_exc);
    title("XZ-Excitation PSF, "  + "Y = 0")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;
    axis image;
    image21.Parent.XLim = [-20,20];
    image21.Parent.YLim = [-20,20];

    subplot(2,3,2);
xzPSF_exc_dither = Data.PSF_exc_3d_dither(:,:,(Data.N+1)/2);
zPSF_exc_dither = squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,(Data.N+1)/2))/max(squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,(Data.N+1)/2)));
image22 = imagesc(X_exc,Z_exc,xzPSF_exc_dither);
     title("Dithered XZ-Excitation PSF, " + ...
          "T_d = " + num2str(Data.dither_period) + "um, " +...
          "Step = " + num2str(Data.dither_step) + ...
          ", Y = 0")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    axis image
    colorbar;
    image22.Parent.XLim = [-20,20];
    image22.Parent.YLim = [-20,20];
    
    subplot(2,3,3);
    hold on
xzOTF_exc = fftshift(fft2(xzPSF_exc));
zOTF_exc = abs(xzOTF_exc(:,(Data.N+1)/2)) / max(abs(xzOTF_exc(:,(Data.N+1)/2)));
% xzOTF_exc_phase = angle(xzOTF_exc(:,(Data.N+1)/2));
image23 = plot( KZ_exc, zOTF_exc);
% phase1 = plot( KZ_exc, xzOTF_exc_phase);
    title("Z-Excitation-OTF, " + "K_X = 0, " + "K_Y = 0")
    ylabel("Normalized a.u. ")
    xlabel("k_z/(4\pin/\lambda_{exc})")
    image23.Color = 'r';
    image23.LineWidth = 2;
    image23.Parent.XLim = [-1,1];
%      phase1.Color = 'g';
%      phase1.LineWidth = 2;
%     lgd = legend("Amplitude","Phase");
    grid on
    hold off

    subplot(2,3,4);
yzPSF_exc = squeeze(Data.PSF_exc_3d(:,(Data.N+1)/2,(Data.N+1)/2:end));
image24 = imagesc(Y_exc((Data.N+1)/2:end), Z_exc, yzPSF_exc);
    title("YZ-Excitation-PSF, " + "X = 0" )
    axis image
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar

    subplot(2,3,5);
yzPSF_exc_dither = squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,(Data.N+1)/2:end));
image25 = imagesc(Y_exc((Data.N+1)/2:end),Z_exc, yzPSF_exc_dither );
    title("Dithered YZ-Excitation-PSF, " + "X = 0" )
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    axis image
    colorbar;

    subplot(2,3,6);
    hold on
    yPSF_exc = squeeze(Data.PSF_exc_3d((Data.N+1)/2,(Data.N+1)/2,:));
    yPSF_exc = yPSF_exc/max(yPSF_exc);
    yPSF_exc_dither = squeeze(Data.PSF_exc_3d_dither((Data.N+1)/2,(Data.N+1)/2,:));
    yPSF_exc_dither = yPSF_exc_dither/max(yPSF_exc_dither);
    % Calculate yFWHM
    index = find(yPSF_exc((Data.N+1)/2:end) <= 0.5);
    index_dither = find(yPSF_exc_dither((Data.N+1)/2:end) <= 0.5);
    if ~isempty(index)
        Data.yFWHM = Y_exc((Data.N+1)/2+index(1)-1);
    else
        Data.yFWHM = "N/A";
    end
    if ~isempty(index_dither)
        Data.yFWHM_dither = Y_exc((Data.N+1)/2+index_dither(1)-1);
    else
        Data.yFWHM = "N/A";
    end
    
image26_1 = plot(Y_exc, yPSF_exc );
image26_2 = plot(Y_exc, yPSF_exc_dither);
    title("Y-Excitation-PSF, " + "X = 0, Z = 0, " + ...
          "yFWHM = " + num2str(Data.yFWHM) + "  \lambda_{exc}/n, " +...
          "yFWHM_{dither} = " + num2str(Data.yFWHM_dither) + "\lambda_{exc}/n")
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("Normalized a.u. ")
    image26_1.Color = 'g';
    image26_1.LineWidth = 2;
    image26_2.Color = 'r';
    image26_2.LineWidth = 2;
    colorbar;
    grid on
    drawnow
    hold off
    
%% Figure 3 - Overall PSF/OTF
    Figures.fig3 = figure(3);  
    Figures.fig3.Name = "Overall Excitation";
    Figures.fig3.WindowState = 'maximized';
    colormap(hot(256))
        
    subplot(2,4,1)
image31 = imagesc(X_exc,Z_exc,Data.PSF_exc_3d_dither(:,:,(Data.N+1)/2) );
    title("Dithered XZ-Excitation PSF")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;
    axis image;
    image31.Parent.XLim = [-10,10];
    image31.Parent.YLim = [-10,10];

    subplot(2,4,2)

 image32 = imagesc(X_det,Z_det,squeeze(Data.PSF_det_3d(:,(Data.N+1)/2,:))');
    title("XZ-Detection PSF ")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;  
    axis image;
    image32.Parent.XLim = [-10,10];
    image32.Parent.YLim = [-10,10];

    subplot(2,4,3)
 overall = squeeze(Data.PSF_det_3d(:,(Data.N+1)/2,:))'.* xzPSF_exc;
 image33 = imagesc(X_exc,Z_exc, overall);
    title("Overall PSF")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;
    axis image;  
    image33.Parent.XLim = [-10,10];
    image33.Parent.YLim = [-10,10];

    subplot(2,4,5)
image35 = imagesc(KX_exc,...
                  KZ_exc,...
                  abs(xzOTF_exc)/max(max(abs(xzOTF_exc))) );
    title("XZ-Excitation OTF ")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image;

    subplot(2,4,6)
 OTF_det_3d = fftshift(fftn(Data.PSF_det_3d));
 xzOTF_det = squeeze(OTF_det_3d(:,(Data.N+1)/2,:))';
 xzOTF_det = abs(xzOTF_det)/max(max(abs(xzOTF_det)));
 image36 = imagesc(KX_det,...
                  KZ_det,...
                  xzOTF_det);
    title("XZ-Detection OTF ")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image;  
    caxis([0,0.5])

    subplot(2,4,7)
xzPSF_det = abs(fftshift(ifft2(squeeze(xzOTF_det))))';
xzPSF_det = xzPSF_det/max(max(xzPSF_det));
Overall_PSF_axial = xzPSF_exc .* xzPSF_det';
Overall_OTF_axial = abs(fftshift(fft2(Overall_PSF_axial)));
Overall_OTF_axial = Overall_OTF_axial/max(max(Overall_OTF_axial));
 image37 = imagesc(KX_det,KZ_det,Overall_OTF_axial);
    title("Overall OTF")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image;  
    caxis([0,0.5])

    h1 = subplot(2,4,[4,8]);
    hold on
line_exc = plot(Z_exc,zPSF_exc_dither);
    line_exc.Color = 'g';
    line_exc.LineWidth = 2;
line_det = plot(Z_det,  squeeze(Data.PSF_det_3d((Data.N+1)/2 ,(Data.N+1)/2,:))');
    line_det.Color = 'b';
    line_det.LineWidth = 2;
line_overall = plot(Z_exc, overall(:,(Data.N+1)/2) );
    line_overall.Color = 'r';
    line_overall.LineWidth = 2;
    title("Overall Axial-PSF")
    xlabel("z/(\lambda_{exc}/n)")
    ylabel("Normalized a.u. ")
    lgd = legend("Excitation", "Detection","Overall");  
        lgd.FontWeight = 'bold';
        lgd.FontSize = 7;
        lgd.LineWidth = 1;
    h1.XLim = [-6,6];
    h1.XTick = linspace(-6,6,13);
    grid on
    hold off
    drawnow
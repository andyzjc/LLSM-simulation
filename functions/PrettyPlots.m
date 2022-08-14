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
    Figures.fig2.Name = "Focal PSF/OTF";
    Figures.fig2.WindowState = 'maximized';
    colormap(hot(256))

     subplot(2,3,1);
xzPSF_exc = Data.PSF_exc_3d(:,:,(Data.N+1)/2);
image21 = imagesc(X_exc,Z_exc,xzPSF_exc);
    title("XZ-Excitation PSF, "  + "Y = 0")
    xlabel("x/  \lambda_{exc}/n")
    ylabel("z/  \lambda_{exc}/n")
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
    xlabel("x/  \lambda_{exc}/n")
    ylabel("z/  \lambda_{exc}/n")
    axis image
    colorbar;
    image22.Parent.XLim = [-20,20];
    image22.Parent.YLim = [-20,20];
    
    subplot(2,3,3);
    hold on
OTF_exc_3d = fftshift(fftn(Data.PSF_exc_3d));
xzOTF_exc = OTF_exc_3d(:,:,(Data.N+1)/2);
zOTF_exc = abs(xzOTF_exc(:,(Data.N+1)/2)) / max(abs(xzOTF_exc(:,(Data.N+1)/2)));
xzOTF_exc_phase = angle(xzOTF_exc(:,(Data.N+1)/2));
image23 = plot( KZ_exc, zOTF_exc);
phase1 = plot( KZ_exc, xzOTF_exc_phase);
    title("Z-Excitation-OTF, " + "K_X = 0, " + "K_Y = 0")
    ylabel("Normalized a.u. ")
    xlabel("kz *   \lambda_{exc}/n")
    image23.Color = 'r';
    image23.LineWidth = 2;
    image23.Parent.XLim = [-1,1];
     phase1.Color = 'g';
     phase1.LineWidth = 2;
    lgd = legend("Amplitude","Phase");
    colorbar;
    grid on
    hold off

    subplot(2,3,4);
yzPSF_exc = squeeze(Data.PSF_exc_3d(:,(Data.N+1)/2,(Data.N+1)/2:end));
image24 = imagesc(Y_exc((Data.N+1)/2:end), Z_exc, yzPSF_exc);
    title("YZ-Excitation-PSF, " + "X = 0" )
    axis image
    xlabel("y/  \lambda_{exc}/n")
    ylabel("z/  \lambda_{exc}/n")
    colorbar

    subplot(2,3,5);
yzPSF_exc_dither = squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,(Data.N+1)/2:end));
image25 = imagesc(Y_exc((Data.N+1)/2:end),Z_exc, yzPSF_exc_dither );
    title("Dithered YZ-Excitation-PSF, " + "X = 0" )
    xlabel("y/  \lambda_{exc}/n")
    ylabel("z/  \lambda_{exc}/n")
    axis image
    colorbar;

    subplot(2,3,6);
    yPSF_exc = squeeze(Data.PSF_exc_3d_dither((Data.N+1)/2,(Data.N+1)/2,:));
    % Calculate yFWHM
    index = find(yPSF_exc >= 0.5);
    Data.yFWHM = Y_exc(index(end)) - Y_exc(index(1)); % pixels
image26 = plot(Y_exc, yPSF_exc );
    title("Y-Excitation-PSF, " + "X = 0, Z = 0, " + ...
          "yFWHM = " + num2str(Data.yFWHM) + "  \lambda_{exc}/n")
    xlabel("y/  \lambda_{exc}/n")
    ylabel("Normalized a.u. ")
    image26.Color = 'r';
    image26.LineWidth = 2;
    colorbar;
    grid on
    drawnow
    
%% Figure 3 - Overall PSF/OTF
    Figures.fig3 = figure(3);  
    Figures.fig3.Name = "Overall XZ-Axial PSF/OTF, focal plane";
    Figures.fig3.WindowState = 'maximized';
    colormap(hot(256))
        
    subplot(2,4,1)
image31 = imagesc(X_exc,Z_exc,Data.PSF_exc_3d_dither(:,:,(Data.N+1)/2) );
    title("Dithered XZ-Excitation PSF")
    xlabel("x/  \lambda_{exc}/n")
    ylabel("z/  \lambda_{exc}/n")
    colorbar;
    axis image;
    image31.Parent.XLim = [-10,10];
    image31.Parent.YLim = [-10,10];

    subplot(2,4,2)
 image32 = imagesc(X_det,Z_det,squeeze(Data.PSF_det_3d(:,(Data.N+1)/2,:))');
    title("XZ-Detection PSF ")
    xlabel("x/  \lambda_{exc}/n")
    ylabel("z/  \lambda_{exc}/n")
    colorbar;  
    axis image;
    image32.Parent.XLim = [-10,10];
    image32.Parent.YLim = [-10,10];

    subplot(2,4,3)
 image33 = imagesc(X_exc,Z_exc,Data.Overall_PSF_axial);
    title("Overall PSF")
    xlabel("x/  \lambda_{exc}/n")
    ylabel("z/  \lambda_{exc}/n")
    colorbar;
    axis image;  
    image33.Parent.XLim = [-10,10];
    image33.Parent.YLim = [-10,10];

    subplot(2,4,5)
image35 = imagesc(KX_exc,...
                  KZ_exc,...
                  abs(xzOTF_exc)/max(max(abs(xzOTF_exc))) );
    title("XZ-Excitation OTF ")
    xlabel("kx *   \lambda_{exc}/n")
    ylabel("kz *   \lambda_{exc}/n")
    colorbar;
    axis image;
    image35.Parent.XLim = [-1,1];
    image35.Parent.YLim = [-1,1];

    subplot(2,4,6)
 OTF_det_3d = fftshift(fftn(Data.PSF_det_3d));
 xzOTF_det = squeeze(OTF_det_3d(:,(Data.N+1)/2,:))';
 xzOTF_det = abs(xzOTF_det)/max(max(abs(xzOTF_det)));
 image36 = imagesc(KX_det,...
                  KZ_det,...
                  xzOTF_det);
    title("XZ-Detection OTF ")
    xlabel("kx *   \lambda_{exc}/n")
    ylabel("kz *   \lambda_{exc}/n")
    colorbar;
    axis image;  
%     image36.Parent.XLim = [-1,1];
%     image36.Parent.YLim = [-1,1];
    caxis([0,0.1])
    subplot(2,4,7)

overallOTF_axial = conv2(abs(xzOTF_exc),abs(xzOTF_det));
overallOTF_axial = abs(overallOTF_axial)/max(max(abs(overallOTF_axial)));
 image37 = imagesc(KX_det,KZ_det,overallOTF_axial);
    title("Overall OTF")
    xlabel("kx *   \lambda_{exc}/n")
    ylabel("kz *   \lambda_{exc}/n")
    colorbar;
    axis image;  
    image37.Parent.XLim = [-1,1];
    image37.Parent.YLim = [-1,1];
    caxis([0,0.1])

    h1 = subplot(2,4,[4,8]);
    hold on
line_exc = plot(zPSF_exc_dither,Z_exc);
    line_exc.Color = 'g';
    line_exc.LineWidth = 2;
line_det = plot(squeeze(Data.PSF_det_3d((Data.N+1)/2,(Data.N+1)/2,:)),Z_det);
    line_det.Color = 'b';
    line_det.LineWidth = 2;
line_overall = plot(squeeze(Data.Overall_PSF_axial(:,(Data.N+1)/2)),Z_exc);
    line_overall.Color = 'r';
    line_overall.LineWidth = 2;
% line_lateral = plot(squeeze(PSF_det_3d(:,(N+1)/2,(N+1)/2)),Z_exc);
%     line_lateral.Color = 'k';
%     line_lateral.LineWidth = 2;
    title("Overall Axial-PSF")
    ylabel("z/  \lambda_{exc}/n")
    xlabel("Normalized a.u. ")
    lgd = legend("Excitation", "Detection","Overall");  
        lgd.FontWeight = 'bold';
        lgd.FontSize = 7;
        lgd.LineWidth = 1;
    h1.YLim = [-6,6];
    h1.YTick = linspace(-6,6,13);
    grid on
    hold off
    drawnow
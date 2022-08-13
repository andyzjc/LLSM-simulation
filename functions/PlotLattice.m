function PlotLattice
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

% Figure 1 - Rear Pupil 
    Figures.fig1 = figure(1);
    Figures.fig1.Name = "XZ-excitation, Y = 0";
    Figures.fig1.WindowState = 'maximized';
    colormap(hot(256))

     subplot(3,4,1);
image11 = imagesc(KX_exc,KZ_exc, real(Lattice.Illumi_ideal) );
    title("Ideal Lattice, " +...
          "N_{xz} = " + num2str(Data.N) + ...
          ", K_{bound} = " + num2str(Data.k_bound) + " (1/um)")
    xlabel("kx *  \lambda/n")
    ylabel("kz *  \lambda/n")
    axis image
    image11.Parent.XLim = [-1,1];
    image11.Parent.YLim = [-1,1];
    colorbar;    

    subplot(3,4,2);
image12 = imagesc(X_exc,Z_exc, abs(Lattice.E_ideal));
    title("Ideal Lattice," + ...
          " \lambda/n_{exc}/n = " + num2str(Data.lambda_exc, '%.3f') + "um / " + ...
           num2str(n))
    xlabel("x/ \lambda/n")
    ylabel("z/ \lambda/n")
    axis image
    image12.Parent.XLim = [-20,20];
    image12.Parent.YLim = [-20,20];
    colorbar;

    subplot(3,4,3);
image13 = imagesc(X_exc, Z_exc, abs(Lattice.E_bound));
    title("Bounded Ideal Lattice, " + ...
          "Bound width = " + num2str(Data.gauss_bound_width) + "um")
    xlabel("x/ \lambda/n")
    ylabel("z/ \lambda/n")
    axis image
    image12.Parent.XLim = [-20,20];
    image12.Parent.YLim = [-20,20];
    colorbar;
    
    subplot(3,4,4)
image14 = imagesc( KX_exc, KZ_exc,...
                  Lattice.Illum_bound );
    title("Rear pupil illumination after bounding")
    xlabel("kx *  \lambda/n")
    ylabel("kz *  \lambda/n")
    axis image
    image14.Parent.XLim = [-1,1];
    image14.Parent.YLim = [-1,1];
    colorbar;

Illum_mask = imfuse(Lattice.Illum_bound,Lattice.A_mask,"falsecolor","ColorChannels","green-magenta");
    subplot(3,4,5)
image15 = imagesc( KX_exc, KZ_exc,...
                  Illum_mask);
    title("Masking, " +...
          "NA_{max} = " + num2str(Data.NAmax) +...
          ", NA_{min} = " + num2str(Data.NAmin) )
    xlabel("kx *  \lambda/n")
    ylabel("kz *  \lambda/n")
    axis image
    image15.Parent.XLim = [-1,1];
    image15.Parent.YLim = [-1,1];

    subplot(3,4,6)
image16 = imagesc( KX_exc, KZ_exc,...
                  Data.Pupil_fun_exc );
    title("Bounded Rear Pupil")
    xlabel("kx *  \lambda/n")
    ylabel("kz *  \lambda/n")
    axis image
    image16.Parent.XLim = [-1,1];
    image16.Parent.YLim = [-1,1];
    colorbar;

    subplot(3,4,7)
xzPSF_exc = Data.PSF_exc_3d(:,:,(Data.N+1)/2);
image17 = imagesc(X_exc, Z_exc, xzPSF_exc);
    title("XZ-Excitation PSF")
    xlabel("x/ \lambda/n")
    ylabel("z/ \lambda/n")
    colorbar;
    axis image;
    image17.Parent.XLim = [-20,20];
    image17.Parent.YLim = [-20,20];

    subplot(3,4,8)
OTF_exc_3d = fftshift(fftn(Data.PSF_exc_3d));
xzOTF_exc = OTF_exc_3d(:,:,(Data.N+1)/2);
image18 = imagesc( KX_exc,...
                  KZ_exc,...
                 abs(xzOTF_exc)) ;
    title("XZ-Excitation OTF")
    xlabel("kx *  \lambda/n")
    ylabel("kz *  \lambda/n")
    colorbar;
    axis image
    image18.Parent.XLim = [-1,1];
    image18.Parent.YLim = [-1,1];

    subplot(3,4,9)
    hold on;
xzPSF_exc_dither = Data.PSF_exc_3d_dither(:,:,(Data.N+1)/2); 
image19 = imagesc(X_exc, Z_exc ,xzPSF_exc_dither);
    title("Dithered XZ-Excitation PSF, " + ...
          "T_d = " + num2str(Data.dither_period) + "um, " +...
          "Dither Step = " + num2str(Data.dither_step))
    xlabel("x/ \lambda/n")
    ylabel("z/ \lambda/n")
    colorbar;
    axis image;
    image19.Parent.XLim = [-10,10];
    image19.Parent.YLim = [-10,10];

    subplot(3,4,10)
zPSF_exc_dither = squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,(Data.N+1)/2))/max(squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,(Data.N+1)/2)));
image110 = plot( Z_exc, zPSF_exc_dither);
    title("Dithered Z-Excitation PSF")
    ylabel("Normalized a.u. ")
    xlabel("z/ \lambda/n")
    image110.LineWidth = 2;
    image110.Color = 'r';
    image110.Parent.YAxis.TickValues = linspace(0,1,11);
    image110.Parent.XLim = [-10,10];
    grid on

    subplot(3,4,11)
OTF_exc_3d_dither = fftshift(fftn(Data.PSF_exc_3d_dither));
xzOTF_exc_dither = OTF_exc_3d_dither(:,:,(Data.N+1)/2);    
image111 = imagesc( KX_exc,...
                    KZ_exc,...
                    abs(xzOTF_exc_dither));
    title("Dithered XZ-Excitation OTF")
    xlabel("kx *  \lambda/n")
    ylabel("kz *  \lambda/n")
    colorbar;
    axis image
    image111.Parent.XLim = [-1,1];
    image111.Parent.YLim = [-1,1];

    subplot(3,4,12)

zOTF_exc_dither = abs(xzOTF_exc_dither(:,(Data.N+1)/2)) / max(abs(xzOTF_exc_dither(:,(Data.N+1)/2)) );
image112 = plot( KZ_exc, zOTF_exc_dither);
    title("Dithered Z-Excitation OTF, " + ...
        "K_x = " + num2str(0) )
    ylabel("Normalized a.u. ")
    xlabel("kz *  \lambda/n")
    image112.LineWidth = 2;
    image112.Color = 'r';
    image112.Parent.XLim = [-1,1];
    image112.Parent.YAxis.TickValues = linspace(0,1,11);
    grid on
    drawnow
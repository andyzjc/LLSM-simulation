function PlotOther
    global Data
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

Figures.fig1 = figure(1);
    Figures.fig1.Name = "XZ-excitation, Y = 0";
    Figures.fig1.WindowState = 'maximized';
    colormap(hot(256))

    subplot(1,3,1)
image16 = imagesc( KX_exc, KZ_exc,...
                  real(Data.Pupil_fun_exc) );
    title("Rear Pupil")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    axis image
    image16.Parent.XLim = [-1,1];
    image16.Parent.YLim = [-1,1];
    colorbar;

    subplot(1,3,2)
xzPSF = Data.PSF_exc_3d(:,:,(N+1)/2);
image17 = imagesc(X_exc, Z_exc, xzPSF);
    title("XZ-Excitation PSF")
    xlabel("x/\lambda")
    ylabel("z/\lambda")
    colorbar;
    axis image;

    subplot(1,3,3)
xzOTF = fftshift(fft2(xzPSF));
image18 = imagesc( KX_exc,...
                  KZ_exc,...
                  abs(xzOTF)/max(max(abs(xzOTF))) ) ;
    title("XZ-Excitation OTF")
    xlabel("kx * \lambda")
    ylabel("kz * \lambda")
    colorbar;
    axis image
    image18.Parent.XLim = [-1,1];
    image18.Parent.YLim = [-1,1];
    drawnow

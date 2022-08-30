function SaveResults
    global Data
    global Lattice
    global Figures

    disp("Saving")
    
    if Lattice.islattice == 1
        foldername = Data.sti_case + "-NAmax=" + num2str(Data.NAmax) + "-NAmin=" + num2str(Data.NAmin);
    else
        foldername = Data.sti_case + "-ApertureNA=" + num2str(Data.apertureNA);
    end

    savingdir = Data.dir + "/" + foldername;
    if ~exist(savingdir, 'dir')
       mkdir(savingdir)
    end
    
    pause(5)
    
    savefig(Figures.fig1, savingdir + "/xzExcitation.fig")
    savefig(Figures.fig2, savingdir + "/yzExcitation.fig")
    savefig(Figures.fig3, savingdir + "/OverallExcitation.fig")

    exportgraphics(Figures.fig1, savingdir + "/xzExcitation.png",'Resolution',500)
    exportgraphics(Figures.fig2, savingdir + "/yzExcitation.png",'Resolution',500)
    exportgraphics(Figures.fig3, savingdir + "/OverallExcitation.png",'Resolution',500)
    
    tic
    for i = (Data.N+1)/2:2:size(Data.Y_exc,2)
        xzPSF_exc = Data.PSF_exc_3d(:,:,i); xzPSF_exc = xzPSF_exc/max(max(xzPSF_exc));
        xzPSF_exc_dither = Data.PSF_exc_3d_dither(:,:,i); xzPSF_exc_dither = xzPSF_exc_dither/max(max(xzPSF_exc_dither));
        zPSF_exc = xzPSF_exc(:,(Data.N+1)/2); zPSF_exc = zPSF_exc/max(zPSF_exc);
        zPSF_exc_dither = xzPSF_exc_dither(:,(Data.N+1)/2); zPSF_exc_dither = zPSF_exc_dither/max(max(zPSF_exc_dither));
        xzOTF_exc = fftshift(fft2(xzPSF_exc)); xzOTF_exc = xzOTF_exc/max(max(xzOTF_exc));
        zOTF_exc = abs(xzOTF_exc(:,(Data.N+1)/2)) / max(abs(xzOTF_exc(:,(Data.N+1)/2)));

        fig4 = figure("Visible","off",'WindowState','maximized');
%         fig4 = figure('WindowState','maximized');
        fig4.Name = "Propagation Profile";
        colormap(hot(256))

        h1 = subplot(4,4,[1:2,5:6]);
        h1_image = imagesc(Data.X_exc,Data.Z_exc,xzPSF_exc);
        title("XZ-Excitation PSF - Y=" + num2str(Data.Y_exc(i),'%.2f') + "\lambda / n")
        h1.Title.FontSize = 5;
        xlabel("x/(\lambda_{exc}/n)")
        ylabel("z/(\lambda_{exc}/n)")
        h1.XAxis.FontSize = 5;
        h1.YAxis.FontSize = 5;
        h1.XAxis.FontWeight = 'bold';
        h1.YAxis.FontWeight = 'bold';
        colorbar;
        axis image;
        h1_image.Parent.XLim = [-20,20];
        h1_image.Parent.YLim = [-20,20];
        

        h2 = subplot(4,4,[3:4,7:8]);
        hold on
        h2_image = plot(Data.Z_exc,zPSF_exc);
        h2_image_dither = plot(Data.Z_exc,zPSF_exc_dither);
        title("Z-Excitation PSF - X=0" + "\lambda / n")
        h2.Title.FontSize = 5;
        ylabel("Normalized a.u. ")
        xlabel("z/(\lambda_{exc}/n)")
        h2.XAxis.FontSize = 4;
        h2.YAxis.FontSize = 5;
        h2.XAxis.FontWeight = 'bold';
        h2.YAxis.FontWeight = 'bold';
        h2_image.LineWidth = 2;
        h2_image.Color = 'r';
        h2_image_dither.LineWidth = 2;
        h2_image_dither.Color = 'g';
        h2_image.Parent.XLim = [-10,10];
        h2_image.Parent.YAxis.TickValues = linspace(0,1,11);
        h2_image.Parent.XAxis.TickValues = linspace(-10,10,21);
        lgd=legend("Non-dither","dithered");
        lgd.FontSize = 5;
        grid on
        axis square
        hold off

        h3 = subplot(4,4,[9:10,13:14]);
        h3_image = imagesc( Data.KX_exc,...
                  Data.KZ_exc,...
                 abs(xzOTF_exc)) ;
        title("XZ-Excitation OTF")
        h3.Title.FontSize = 5;
        xlabel("k_x/(4\pin/\lambda_{exc})")
        ylabel("k_z/(4\pin/\lambda_{exc})")
        h3.XAxis.FontSize = 5;
        h3.YAxis.FontSize = 5;
        h3.XAxis.FontWeight = 'bold';
        h3.YAxis.FontWeight = 'bold';
        colorbar;
        axis image
        h3_image.Parent.XLim = [-0.5,0.5];
        h3_image.Parent.YLim = [-0.5,0.5];

        h4 = subplot(4,4,[11:12,15:16]);
        h4_image = plot(Data.KZ_exc,zOTF_exc);
        title("Z-Excitation OTF, K_X=0")
        h4.Title.FontSize = 5;
        xlabel("kz * \lambda / n")
        ylabel("Normalized a.u. ")
        h4.XAxis.FontSize = 5; 
        h4.YAxis.FontSize = 5; 
        h4.XAxis.FontWeight = 'bold'; 
        h4.YAxis.FontWeight = 'bold'; 
        h4_image.Color = 'r';
        h4_image.LineWidth = 2;
        h4_image.Parent.XLim = [-1,1];
        h4_image.Parent.YAxis.TickValues = linspace(0,1,11);
        h4.XLim = [-0.5,0.5];
        h4_image.Parent.XAxis.TickValues = linspace(-0.5,0.5,11);
        grid on
        axis square

        exportgraphics(fig4,savingdir + "/profile.gif",'Append',true,'Resolution',500)
        clear fig4
        delete fig4
    end
    toc
    
    disp("Saved to " + savingdir)
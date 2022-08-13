function SaveResults
    global Data
    global Figures

    disp("Saving")
    
    foldername = [];

    savingdir = Data.dir + "/" + foldername;
    if ~exist(savingdir, 'dir')
       mkdir(savingdir)
    end
    
    save(savingdir + '/results.mat','-struct','Data','-v7.3');
    savefig(Figures.fig1, savingdir + '/XZ-Excitation.fig')
    savefig(Figures.fig2, savingdir + '/PSF-OTF.fig')
    savefig(Figures.fig3,savingdir + '/Overall.fig')

    saveas(Figures.fig1, savingdir + '/XZ-Excitation.png')
    saveas(Figures.fig2, savingdir + '/PSF-OTF.png')
    saveas(Figures.fig3, savingdir + '/Overall.png')
    
    tic

    for i = (Data.N+1)/2:size(Data.Y_exc,2)
        Figures.fig4 = figure("Visible","off",'WindowState','maximized');
        Figures.fig4.Name = "PSF/OTF";
        h1 = subplot(2,2,1);
        imagesc(Data.X_exc,Data.Z_exc,squeeze(Data.PSF_exc_3d(:,:,i)));
        title("XZ-Excitation PSF - Y=" + num2str(Data.Y_det(i)) + "\lambda / n")
        xlabel("x/\lambda / n")
        ylabel("z/\lambda / n")
        h1.XAxis.FontSize = 15;
        h1.YAxis.FontSize = 15;
        h1.XAxis.FontWeight = 'bold';
        h1.YAxis.FontWeight = 'bold';
        colormap(hot(256))
        colorbar;
        axis image;

        h2 = subplot(2,2,2);
        imagesc(Data.KX_exc, Data.KZ_exc, squeeze(Data.OTF_exc_3d(:,:,i)));
        title("XZ-Excitation OTF - Y=" + num2str(Data.Y_det(i)) + "\lambda / n")
        xlabel("kx * \lambda / n")
        ylabel("kz * \lambda / n")
        h2.XAxis.FontSize = 15; 
        h2.YAxis.FontSize = 15; 
        h2.XAxis.FontWeight = 'bold'; 
        h2.YAxis.FontWeight = 'bold'; 
        colormap(hot(256))
        colorbar; 
        axis image; 

        h3 = subplot(2,2,3);
        if Data.Dither == 1
            line = plot(Data.Z_exc, squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,i)));
            title("Dithered Z-Excitation PSF - Y=" + num2str(Data.Y_det(i)) + "\lambda / n")
        else
            line = plot(Data.Z_exc, squeeze(Data.PSF_exc_3d(:,(Data.N+1)/2,i)));
            title("Z-Excitation PSF - Y=" + num2str(Data.Y_det(i)) + "\lambda / n")
        end
        xlabel("z/\lambda / n")
        ylabel("Normalized a.u. ")
        h3.XAxis.FontSize = 15; 
        h3.YAxis.FontSize = 15; 
        h3.XAxis.FontWeight = 'bold'; 
        h3.YAxis.FontWeight = 'bold'; 
        line.Color = 'r';
        line.LineWidth = 2;
        colorbar; 
        grid on;

        h4 = subplot(2,2,4);
        hold on
        if Data.Dither == 1
            amp = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d_dither(:,(Data.N+1)/2,i)));
            phase = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d_dither_phase(:,(Data.N+1)/2,i)));
            title("Dithered Z-Excitation OTF - Y=" + num2str(Data.Y_det(i)) + "\lambda / n")
        else
            amp = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d(:,(Data.N+1)/2,i)));
            phase = plot(Data.KZ_exc, squeeze(Data.OTF_exc_3d_phase(:,(Data.N+1)/2,i)));
            title("Z-Excitation OTF - Y=" + num2str(Data.Y_det(i)) + "\lambda / n")
        end
        xlabel("kz * \lambda / n")
        ylabel("Normalized a.u. ")
        h4.XAxis.FontSize = 15; 
        h4.YAxis.FontSize = 15; 
        h4.XAxis.FontWeight = 'bold'; 
        h4.YAxis.FontWeight = 'bold'; 
        amp.Color = 'r';
        amp.LineWidth = 2;
        phase.Color = 'g';
        phase.LineWidth = 2;
        colorbar; 
        hold off
        grid on;

        frame = getframe(Figures.fig4);
        imwrite(frame2im(frame),hot, savingdir + '/XZ-PSF-OTF.tif','WriteMode','append')
    end
    toc

    disp("Saved to " + savingdir)
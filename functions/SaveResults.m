function SaveResults
    global Data
    global Lattice
    global Figures

    disp("Saving")
    
    if Lattice.islattice == 1
        foldername = Data.sti_case + ",NAmax=" + num2str(Data.NAmax) + ",NAmin=" + num2str(Data.NAmin);
    else
        foldername = Data.sti_case + ",ApertureNA=" + num2str(Data.apertureNA);
    end

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
    OTF_exc_3d = fftshift(fftn(Data.PSF_exc_3d));
    OTF_exc_3d_dither = fftshift(fftn(Data.PSF_exc_3d_dither));
    for i = (Data.N+1)/2:size(Data.Y_exc,2)
        Figures.fig4 = figure("Visible","off",'WindowState','maximized');
        Figures.fig4.Name = "PSF/OTF";
        h1 = subplot(2,3,1);
        imagesc(Data.X_exc,Data.Z_exc,squeeze(Data.PSF_exc_3d(:,:,i)));
        title("XZ-Excitation PSF - Y=" + num2str(Data.Y_det(i),'%.2f') + "\lambda / n")
        xlabel("x/\lambda / n")
        ylabel("z/\lambda / n")
        h1.XAxis.FontSize = 10;
        h1.YAxis.FontSize = 10;
        h1.XAxis.FontWeight = 'bold';
        h1.YAxis.FontWeight = 'bold';
        colormap(hot(256))
        colorbar;
        axis image;

        h2 = subplot(2,3,2);
        image1 = imagesc(Data.KX_exc, Data.KZ_exc, squeeze(abs(OTF_exc_3d(:,:,i))));
        title("XZ-Excitation OTF - Y=" + num2str(Data.Y_det(i),'%.2f') + "\lambda / n")
        xlabel("kx * \lambda / n")
        ylabel("kz * \lambda / n")
        h2.XAxis.FontSize = 10; 
        h2.YAxis.FontSize = 10; 
        h2.XAxis.FontWeight = 'bold'; 
        h2.YAxis.FontWeight = 'bold'; 
        axis image;
        image1.Parent.XLim = [-1,1];
        image1.Parent.YLim = [-1,1];        
        colormap(hot(256))
        colorbar; 
         

        h3 = subplot(2,3,3);
        line = plot(Data.Z_exc, squeeze(Data.PSF_exc_3d_dither(:,(Data.N+1)/2,i)));
        title("Dithered Z-Excitation PSF - Y=" + num2str(Data.Y_det(i),'%.2f') + "\lambda / n")
        xlabel("z/\lambda / n")
        ylabel("Normalized a.u. ")
        h3.XAxis.FontSize = 10; 
        h3.YAxis.FontSize = 10; 
        h3.XAxis.FontWeight = 'bold'; 
        h3.YAxis.FontWeight = 'bold'; 
        line.Color = 'r';
        line.LineWidth = 2;
        colorbar; 
        grid on;

        h4 = subplot(2,3,4);
        amp = plot(Data.KZ_exc, squeeze(abs(OTF_exc_3d_dither(:,(Data.N+1)/2,i)/max(abs(OTF_exc_3d_dither(:,(Data.N+1)/2,i))))));
        title("Dithered Z-Excitation OTF - Y=" + num2str(Data.Y_det(i),'%.2f') + "\lambda / n")
        xlabel("kz * \lambda / n")
        ylabel("Normalized a.u. ")
        h4.XAxis.FontSize = 10; 
        h4.YAxis.FontSize = 10; 
        h4.XAxis.FontWeight = 'bold'; 
        h4.YAxis.FontWeight = 'bold'; 
        amp.Color = 'r';
        amp.LineWidth = 2;
        grid on;

        h5 = subplot(2,3,5:6);
        phase = plot(Data.KZ_exc,  squeeze(angle(OTF_exc_3d_dither(:,(Data.N+1)/2,i))));
        title("Dithered Z-Excitation OTF Phase - Y=" + num2str(Data.Y_det(i),'%.2f') + "\lambda / n")
        xlabel("kz * \lambda / n")
        ylabel("Normalized a.u. ")
        h4.XAxis.FontSize = 10; 
        h4.YAxis.FontSize = 10; 
        h4.XAxis.FontWeight = 'bold'; 
        h4.YAxis.FontWeight = 'bold'; 
        phase.Color = 'r';
        phase.LineWidth = 2;
        phase.Parent.YLim = [-pi,pi];
        grid on;

        frame = getframe(Figures.fig4);
        imwrite(frame2im(frame),hot, savingdir + '/XZ-PSF-OTF.tif','WriteMode','append')
    end
    toc

    disp("Saved to " + savingdir)
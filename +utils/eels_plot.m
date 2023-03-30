function eels_plot(e_w, t_w, psi)

    imagesc(e_w,t_w, psi);
    ylim([-1 , 1.5]);
    xlim([-5,5])
    colormap jet
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.YTick = -1:0.5:1.5;
    ax.XTick = -4:2:4;
    ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    

end

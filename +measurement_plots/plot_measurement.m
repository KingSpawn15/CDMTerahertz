function plot_measurement(saved_measurement_file_name, measurement_polarization)

load(saved_measurement_file_name);

% find index closest to polarization angle
[~, ind] = min( abs( Polarization - measurement_polarization ) );

RepetitionInd = 1;


for PolInd = ind
    
    EnergyCrop = squeeze(EnergyCropAll(RepetitionInd,PolInd,:))';
    DataSetCrop = squeeze(DataSetCropAll(RepetitionInd,PolInd,:,:))';
    
    %         dataset = movmean(DataSetCrop(:,2:end)',3,1);
    dataset = DataSetCrop(:,2:end);
    dataset = dataset./trapz(EnergyCrop,dataset,1);
    
    figure(7)
    imagesc(EnergyCrop,time,dataset')
    xlabel('Energy shift [eV]')
    ylabel('Time delay [ps]')
    xlim([min(EnergyCrop) max(EnergyCrop)])
    ylim([min(time) 1.5])
    colormap jet
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.YTick = -1:0.5:2.5;
    box on
    caxis([0 max(dataset(:))])
    colorbar
    annotation('textbox',[0.2, 0.2, 0.1, 0.1],'String',['\theta = ',num2str(Polarization(PolInd)),'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',18)
    drawnow
  
end

end



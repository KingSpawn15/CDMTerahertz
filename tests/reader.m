
% aa = squeeze(cell2mat(DataSetCropAll(11 )));
deltat = Time;
energy = EnergyCrop;
% imagesc(energy, deltat, aa')

close all

tiledlayout(1,6);

for i = 1:2:11
    
    if i ~= 1 || i~=11
        nexttile
    end
   
    eels = squeeze(cell2mat(DataSetCropAll(i)));
    imagesc(energy, deltat, eels')
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
% 
% close all
% 
% tiledlayout(1,11);
% 
% for i = 1:11
%     
%     
% 
%     
% 
% 
% end
% % 
% plot(1:10,1:10)
% 
% % Tile 2
% nexttile
% plot(1:10,1:10)
% 
% % Tile 3
% nexttile
% plot(1:10,1:10)
% 
% % Tile 4
% nexttile
% plot(1:10,1:10)
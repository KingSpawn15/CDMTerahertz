% aa = squeeze(cell2mat(DataSetCropAll(11 )));
% the above line is commented out, so it is not used in this code
load('saved_matrices\PulseEnergy.mat')
deltat = Time;
energy = EnergyCrop;
% assign variables deltat and energy to values of Time and EnergyCrop, respectively
% these variables will be used later in the code

close all
% close all previously opened figures

tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'tight');

% create a 1 row and 6 column tiled layout for the plots that follow

for i = [1,5,10]
% loop over the odd-numbered indices from 1 to 11 (i.e., 1, 3, 5, 7, 9, 11)

    if i ~= 1 || i~=11
    % if the current index i is not 1 or 11

        nexttile
        % select the next tile in the layout (i.e., move to the next subplot)

    end
   
    eels = squeeze(cell2mat(DataSetCropAll(i)));

    imagesc(energy, deltat, eels')

    ylim([-.5 , 1.5]);
    xlim([-4.5,4.5])
    % set limits for the y-axis and x-axis, respectively

    colormap jet
    % set the colormap to 'jet'

    axis square
    % set the aspect ratio of the plot to 1:1

    ax = gca;
    ax.FontSize = 20;
    ax.FontName = 'helvetica'
    ax.LineWidth = 1;
    
    ax.YTick = []
    ax.XTick = -4:2:4;
    % modify the current axis properties, including fontsize, linewidth,
    % y-axis ticks, and x-axis ticks
    if i == 1
    ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',22, 'FontName' , 'helvetica');
%     ax.YTick = -1:0.5:1.5;
    ax.YTick = -.5:0.5:1.5;
    end
    xlabel('Energy [eV]','Color',[0.3 0.3 0.3],'FontSize',22, 'FontName' , 'helvetica');
    % add y-label and x-label to the plot, with specified color and fontsize
    
end
% set(gcf,'Renderer','painters')
set(gcf, 'Position', get(gcf, 'Position') .* [1 1 2 1]);
% set the position of the figure to be wider than its height, with a width-to-height ratio of 5:1

% Export figure as SVG
% print(gcf, 'results/power_90degree.svg');

% Export figure as PNG with resolution of 300 DPI
exportgraphics(gcf, 'cleo_plots/aux_figures/power_90degree.png', 'Resolution', 300, 'BackgroundColor', 'none');
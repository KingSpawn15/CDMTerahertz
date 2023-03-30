% aa = squeeze(cell2mat(DataSetCropAll(11 )));
% the above line is commented out, so it is not used in this code
deltat = Time;
energy = EnergyCrop;
% assign variables deltat and energy to values of Time and EnergyCrop, respectively
% these variables will be used later in the code

close all
% close all previously opened figures

tiledlayout(1,6);
% create a 1 row and 6 column tiled layout for the plots that follow

for i = 1:2:11
% loop over the odd-numbered indices from 1 to 11 (i.e., 1, 3, 5, 7, 9, 11)

    if i ~= 1 || i~=11
    % if the current index i is not 1 or 11

        nexttile
        % select the next tile in the layout (i.e., move to the next subplot)

    end
   
    eels = squeeze(cell2mat(DataSetCropAll(i)));
    % assign variable eels to the data from DataSetCropAll at index i
    % the data is squeezed and converted to a matrix using cell2mat

    imagesc(energy, deltat, eels')
    % create an image plot of the data
    % energy is used for the x-axis, deltat for the y-axis, and eels for the data
    % the transpose of eels is used to orient the plot correctly

    ylim([-1 , 1.5]);
    xlim([-5,5])
    % set limits for the y-axis and x-axis, respectively

    colormap jet
    % set the colormap to 'jet'

    axis square
    % set the aspect ratio of the plot to 1:1

    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.YTick = -1:0.5:1.5;
    ax.XTick = -4:2:4;
    % modify the current axis properties, including fontsize, linewidth,
    % y-axis ticks, and x-axis ticks

    ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    % add y-label and x-label to the plot, with specified color and fontsize
    
end
% set(gcf,'Renderer','painters')
set(gcf, 'Position', get(gcf, 'Position') .* [1 1 5 1]);
% set the position of the figure to be wider than its height, with a width-to-height ratio of 5:1

% Export figure as SVG
% print(gcf, 'results/power_90degree.svg');

% Export figure as PNG with resolution of 300 DPI
exportgraphics(gcf, 'results/power_90degree.png', 'Resolution', 300, 'BackgroundColor', 'none');
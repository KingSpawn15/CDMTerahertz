% aa = squeeze(cell2mat(DataSetCropAll(11 )));
% the above line is commented out, so it is not used in this code
load('saved_matrices\PulseEnergy.mat')
deltat = Time;
energy = EnergyCrop;
% assign variables deltat and energy to values of Time and EnergyCrop, respectively
% these variables will be used later in the code

close all

x_c = zeros(11, size(Time,2));

for ii = 1:11
    eels = squeeze(cell2mat(DataSetCropAll(ii)))';
    eels = eels ./ sqrt(sum(eels.^2, 2));
    % Define the mass of each element in the matrix (assuming unit mass)
    % mass_matrix = ones(size(eels));
    mass_matrix = eels;
    % Initialize a vector to store the x coordinates of the center of mass of each row
   

    % Calculate the x coordinate of the center of mass of each row
    for i = 1:size(eels,1)
        row_mass = sum(mass_matrix(i,:));
        for j = 1:size(eels,2)
            x_c(ii,i) = x_c(ii,i) + j*mass_matrix(i,j);
        end
        x_c(ii,i) = x_c(ii,i)/row_mass;
    end
end

com = energy(floor(x_c));

% close all previously opened figures

tiledlayout(1,8);
% create a 1 row and 6 column tiled layout for the plots that follow

for i = [1:2:9,10,11]
% loop over the odd-numbered indices from 1 to 11 (i.e., 1, 3, 5, 7, 9, 11)

    if i ~= 1 || i~=11
    % if the current index i is not 1 or 11

        nexttile
        % select the next tile in the layout (i.e., move to the next subplot)

    end
   
    eels = squeeze(cell2mat(DataSetCropAll(i)))';
    % assign variable eels to the data from DataSetCropAll at index i
    % the data is squeezed and converted to a matrix using cell2mat

    imagesc(energy, deltat, eels)
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

%%
% close all
% figure;

colors = turbo(10);
for i = 1:1:10
    plot(com(i,:) + 2 * (i-1), deltat, 'color', colors(i,:),'LineWidth',1);
%     plot( deltat .* 0 + 2 * (i-1), deltat, 'color', colors(i,:),'LineStyle','--', ...
%         'LineWidth',.5)
set(gca, 'YDir', 'reverse');
    hold on;
end
% axis equal;
% daspect([1 1 1]);
ylim([-.5 , 2]);
xlim([-2,22])
ax = gca;
ax.FontSize = 14;
% ax.LineWidth = 2;
ax.YTick = -1:0.5:2;
ax.XTick = [];

stretch_factor = 1.5; % the desired stretching factor
pos = get(gcf, 'Position'); % get the current position vector of the figure
pos(2) = pos(2)*0.5
pos(3) = pos(3) * stretch_factor; % stretch the height by the stretching factor
pos(4) = pos(4) / stretch_factor;
set(gcf, 'Position', pos); % set the new position vector of the figure
ylabel('\Delta t [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
% Export figure as PNG with resolution of 300 DPI
exportgraphics(gcf, 'cleo_plots/aux_figures/power_90degree.png', 'Resolution', 300, 'BackgroundColor', 'none');
% 
% eels = eels ./ sqrt(sum(eels.^2, 2));
% % Define the mass of each element in the matrix (assuming unit mass)
% % mass_matrix = ones(size(eels));
% mass_matrix = eels;
% % Initialize a vector to store the x coordinates of the center of mass of each row
% x_c = zeros(1, size(eels,1));
% 
% % Calculate the x coordinate of the center of mass of each row
% for i = 1:size(eels,1)
%     row_mass = sum(mass_matrix(i,:));
%     for j = 1:size(eels,2)
%         x_c(i) = x_c(i) + j*mass_matrix(i,j);
%     end
%     x_c(i) = x_c(i)/row_mass;
% end
% 
% 
% plot(deltat, energy(floor(x_c)))
% xlim([-1 , 1.5]);
% ylim([-5,5])
% ax = gca;
% ax.FontSize = 18;
% ax.LineWidth = 1;
% ax.XTick = -1:0.5:1.5;
% ax.YTick = -4:2:4;

set(gcf,'Renderer','painters')
% set(gcf, 'Position', get(gcf, 'Position') .* [1 1 7 1]);
% % set the position of the figure to be wider than its height, with a width-to-height ratio of 5:1
% 
% % Export figure as SVG
% print(gcf, 'results/power_90degree.svg');
% 
% % Export figure as PNG with resolution of 300 DPI
exportgraphics(gcf, 'results/power_90degree.png', 'Resolution', 300, 'BackgroundColor', 'none');
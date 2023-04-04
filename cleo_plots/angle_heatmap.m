theta_chart = [];
for angle = 0 : 1: 178
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    offset = 100;
    theta_chart = vertcat(theta_chart,[mod(fix(2*(-offset + angle)),360),angle]);

end
theta_chart = sortrows(theta_chart);

data = [];

for angle = [0: 2: 178]
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
    disp(angle)
    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);
    data = cat(3,data,dataset);

end

x_c = zeros(90, size(time,2));
for ii = 1:90
    eels = data(:,:,ii);
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

%%

close all

tiledlayout(3,1)
energy = EnergyCrop
com = energy(floor(x_c));

ax1 = nexttile
imagesc(time, 0:2:178, com, [-3, 5])

% Set the font size for the axis labels and title
set(gca, 'FontSize', 10);
colormap(ax1,'jet');
title('My Title', 'FontSize', 12);
xlabel('X Axis Label', 'FontSize', 10);
ylabel('Y Axis Label', 'FontSize', 10);
axis  square;

ax2 = nexttile
imagesc(time, 0:2:178, com, [-3, 5])
colormap(ax2,'hot');
% Set the font size for the axis labels and title
set(gca, 'FontSize', 10);
title('My Title', 'FontSize', 12);
xlabel('X Axis Label', 'FontSize', 10);
ylabel('Y Axis Label', 'FontSize', 10);
axis  square;

ax3 = nexttile
imagesc(time, 0:2:178, com, [-3, 5])
colormap(ax3,'copper');
% Set the font size for the axis labels and title
set(gca, 'FontSize', 10);
title('My Title', 'FontSize', 12);
xlabel('X Axis Label', 'FontSize', 10);
ylabel('Y Axis Label', 'FontSize', 10);
axis  square;

pos = get(gcf, 'Position')
set(gcf, 'Position', pos .* [1 .1 1 2])

% ii = 1;
% for angle = [0: 1: 178]
%     %     angle = 10;
%     %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
% 
%     [dataset, EnergyCrop,...
%         time] = measurement_plots.data_measurement(angle);
% 
%     offset = 100;
%     theta = mod(fix(2*(-offset + angle)),360);
% 
%     if ~((mod(theta,10) == 0) || (mod(theta+2,10) == 0))
%         continue
%     end
% 
%     theta_str = num2str(theta);
% 
%     tc = time(time > -1 & time < 1.5);
%     psi = dataset(time > -1 & time < 1.5,:);
% 
%     ii = ii + 1;
%%
%     close all;
%     figure;
% 
%     imagesc(EnergyCrop,tc,psi);
%     annotation('textbox',[0.15, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
%         'color',[1 1 1],'LineStyle','none','FontSize',18);
%     annotation('textbox',[0.15, 0.1, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
%         'color',[1 1 1],'LineStyle','none','FontSize',8);
%     colormap('bone');
%     colorbar;
%     set(gca,'FontSize',16);
%     xlabel('Energy shift [eV]','FontSize',18,'FontName','Times New Roman');
%     ylabel('Time delay [ps]','FontSize',18,'FontName','Times New Roman');
%     yticks([-1 : 0.2 : 1.5]);
% 
%     xticks([-4:2:4]);
%     ylim([-1,1.5]);
%     xlim([-5,5]);
% 
% %%
%     str = [dir,...
%         'single_theta=',num2str(theta),...
%         ];
%     exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
% 
%     disp([angle,theta]);
% end

clear all;
close all;
dir = 'examples/results/cleo_plots/';
[~,~,~] = mkdir(dir);
ii = 1;


theta_chart = [];
for angle = 0 : 1: 178
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    offset = 102;
    theta_chart = vertcat(theta_chart,[mod(fix(2*(-offset + angle)),360),angle]);

end
theta_chart = sortrows(theta_chart);

ii = 1;
for angle = [0: 1: 178]
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);

    offset = 102;
    theta = mod(fix(2*(-offset + angle)),360);

    if ~((mod(theta,10) == 0) || (mod(theta+2,10) == 0))
        continue
    end

    theta_str = num2str(theta);

    tc = time(time > -1 & time < 1.5);
    psi = dataset(time > -1 & time < 1.5,:);

    ii = ii + 1;
%%
    close all;
    figure;

    imagesc(EnergyCrop,tc,psi);
    annotation('textbox',[0.15, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',18);
    annotation('textbox',[0.15, 0.1, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',8);
    colormap('bone');
    colorbar;
    set(gca,'FontSize',16);
    xlabel('Energy shift [eV]','FontSize',18,'FontName','Times New Roman');
    ylabel('Time delay [ps]','FontSize',18,'FontName','Times New Roman');
    yticks([-1 : 0.2 : 1.5]);

    xticks([-4:2:4]);
    ylim([-1,1.5]);
    xlim([-5,5]);

%%
    str = [dir,...
        'single_theta=',num2str(theta),...
        ];
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);

    disp([angle,theta]);
end

%%

for angle = [107 147]
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);

    offset = 102;
    theta = mod(fix(2*(-offset + angle)),360);

    if ~((mod(theta,10) == 0) || (mod(theta+2,10) == 0))
        continue
    end

    theta_str = num2str(theta);

    tc = time(time > -1 & time < 1.5);
    psi = dataset(time > -1 & time < 1.5,:);

    ii = ii + 1;
%%
    close all;
    figure;
    
    imagesc(EnergyCrop,tc,psi);
  ax1 = gca;
     axis square;
%     annotation('textbox',[0.15, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
%         'color',[1 1 1],'LineStyle','none','FontSize',18);
%     annotation('textbox',[0.15, 0.1, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
%         'color',[1 1 1],'LineStyle','none','FontSize',8);
    colormap('hot');
%     colorbar;
%     set(gca,'FontSize',16);
%     xlabel('Energy shift [eV]','FontSize',18,'FontName','Times New Roman');
%     ylabel('Time delay [ps]','FontSize',18,'FontName','Times New Roman');
set(ax1.XAxis,'FontName' , 'arial','FontSize',16);
    set(ax1.YAxis,'FontName' , 'arial','FontSize',16);
    xlabel('Energy shift (eV)','FontName','arial');
    ylabel('\Delta t delay (ps)','FontName','arial');
    yticks([-1 : 0.5 : 1.5]);

    xticks([-4:2:4]);
    ylim([-1,1.5]);
    xlim([-5,5]);

%%
    str = [dir,...
        'single_theta=',num2str(theta),...
        ];
    set(gcf,'Position',[100, 50, 300, 300]);
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);

    disp([angle,theta]);
end




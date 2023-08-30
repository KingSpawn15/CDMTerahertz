clear all;
close all;
dir = 'results/fourier_analysis/measurement_plots/';
[~,~,~] = mkdir(dir);
ii = 1;

offset = 106;
theta_chart = [];
for angle = 0 : 1: 178
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    
    theta_chart = vertcat(theta_chart,[mod(fix(2*(-offset + angle)),360),angle]);

end
theta_chart = sortrows(theta_chart);

angle_p2 =  theta_chart((theta_chart(:,1)==44 |  theta_chart(:,1)==134),2).';

angle_p =  theta_chart((theta_chart(:,1)==44 |  theta_chart(:,1)==134) |...
    (mod(theta_chart(:,1),10) == 0 &...
    theta_chart(:,1) <= 180&...
    theta_chart(:,1) >=10),2).';
close all;

nrows = 2;
ncol = 10;
figure;
tiledlayout(2,10,'TileSpacing','compact');
ii = 0;
for angle = angle_p
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);

%     offset = 100;
    theta = mod(fix(2*(-offset + angle)),360);

%     if ~((mod(theta,10) == 0) || (mod(theta+2,10) == 0))
%         continue
%     end

    theta_str = num2str(theta);

    tc = time(time > -1 & time < 1.5);
    psi = dataset(time > -1 & time < 1.5,:);

    ii = ii + 1;
%%
    
    nexttile;

    imagesc(EnergyCrop,tc,psi);
%     annotation('textbox',[0.15, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
%         'color',[1 1 1],'LineStyle','none','FontSize',18);
%     annotation('textbox',[0.15, 0.1, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
%         'color',[1 1 1],'LineStyle','none','FontSize',8);
    colormap('bone');
%     colorbar;
%     set(gca,'FontSize',16);
    if ii > ncol
%         xlabel('Energy shift [eV]','FontSize',18,'FontName','Times New Roman');
    end

    if (ii == 1 || ii ==ncol + 1)
%     ylabel('Time delay [ps]','FontSize',18,'FontName','Times New Roman');
    end

        yticks([]);

    xticks([]);

%     yticks([-1 : 0.2 : 1.5]);

%     xticks([-4:2:4]);
    ylim([-1,1.5]);
    xlim([-5,5]);

%%
%     str = [dir,...
%         'single_theta=',num2str(theta),...
%         ];
%     exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);

    disp([angle,theta]);
end

str = [dir,...
        'combined_offset_102',...
        ];
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
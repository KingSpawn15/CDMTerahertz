clear all;
close all;
dir = 'examples/results/cleo_plots/';
[~,~,~] = mkdir(dir);
ii = 1;

offset = 102;
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

% angle_p =  theta_chart(((theta_chart(:,1)==10 |  theta_chart(:,1)==30) |...
%     (theta_chart(:,1)==44 |  theta_chart(:,1)==70) |...
%     (theta_chart(:,1)==90 |  theta_chart(:,1)==100) |...
%     (theta_chart(:,1)==120 |  theta_chart(:,1)==134) |...
%     (theta_chart(:,1)==160 |  theta_chart(:,1)==180) &...
%     theta_chart(:,1) <= 180 &...
%     theta_chart(:,1) >=10),2).';

angle_p =  theta_chart(((theta_chart(:,1)==10 |  theta_chart(:,1)==30) |...
    (theta_chart(:,1)==50 |  theta_chart(:,1)==70) |...
    (theta_chart(:,1)==90 |  theta_chart(:,1)==110) |...
    (theta_chart(:,1)==130 |  theta_chart(:,1)==150) |...
    (theta_chart(:,1)==170)  &...
    theta_chart(:,1) <= 180 &...
    theta_chart(:,1) >=10),2).';

close all;

nrows = 1;
ncol = 5;
figure;
tiledlayout(1,10,'TileSpacing','compact','Padding','compact');
ii = 0;
for angle = angle_p

    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);


    theta = mod(fix(2*(-offset + angle)),360);

    theta_str = num2str(theta);

    tc = time(time > -1 & time < 1.5);
    psi = dataset(time > -1 & time < 1.5,:);

    ii = ii + 1;
    ax1 = nexttile;
    
    imagesc(EnergyCrop,tc,psi);
   pbaspect(ax1,[1 1 1]); % <---- move to after-plot
    %     annotation('textbox',[0.15, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
    %         'color',[1 1 1],'LineStyle','none','FontSize',18);
    %     annotation('textbox',[0.15, 0.1, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
    %         'color',[1 1 1],'LineStyle','none','FontSize',8);
    colormap('hot');
    set(ax1.XAxis,'FontName' , 'arial','FontSize',12);
     set(ax1.YAxis,'FontName' , 'arial','FontSize',12);
   
    if ii == 1
        ylabel('\Delta t delay (ps)','FontName','arial');
         
    end

    if ii == 5
       
         xlabel('Energy shift (eV)','FontName','arial');
    end

    %     colorbar;
    %     set(gca,'FontSize',16);
    yticks([]);
    xticks([]);
    if ii > ncol
%         xlabel('Energy shift [eV]','FontName','Times New Roman');
        xticks([-4:2:4]);
    end
xticks([-4:2:4]);
    if (ii == 1 )
%         ylabel('\Delta t (ps)','FontName','Times New Roman');
        yticks([-1 : 0.5 : 1.5]);
    end

    ylim([-1,1.5]);
    xlim([-5,5]);

    %%
    %     str = [dir,...
    %         'single_theta=',num2str(theta),...
    %         ];
    %     exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);

    disp([angle,theta]);
end
%     axis equal;            % <---- move to after-plot
%     daspect(ax1,[1 1 1]);  % <---- move to after-plot
 
str = [dir,...
    'combined_cleo_measurement_straight',...
    ];
set(gcf,'Position',[50,100,1600,200]);
exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);








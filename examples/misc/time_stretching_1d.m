clear all;
close all;
dir = 'results/stretching_square/oned/';
[~,~,~] = mkdir(dir);
ii = 1;


for angle = 0 : 2 : 178
    %     angle = 10;
    % [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
    
    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);
    
    offset = 100;
    theta = mod(fix(2*(-offset + angle)),360);
    theta_str = num2str(theta);
    
    tc = time(time < -0.5 | time > 0.8);
    dc = dataset(time < -0.5 | time > 0.8,:);
    [me, ind] = max(dc,[],2) ;
    eind = EnergyCrop(ind);
    
    tc_all = time(time > -1 & time < 1.5);
    dc_all = dataset(time > -1 & time < 1.5,:);
    [me_all, ind_all] = max(dc_all,[],2) ;
    eind_all = EnergyCrop(ind_all);
    
    
    [fitresult, gof] = fit(tc(:), eind(:), 'poly1');
    fitdata = fitresult(tc_all);
    
    
    plot(eind_all, tc_all ,'LineStyle','none',...
        'Marker',"o",...
        'MarkerFaceColor',[0.8,0,0],...
        'MarkerSize',5,...
        'Color',[0.8,0,0]);
    hold on;
    plot(fitdata, tc_all,'k--');
    set(gca,'FontSize',12);
    ylabel('\Delta t (ps)','FontSize', 14);
    xlabel('Energy (ev)','FontSize', 14);
    xlim([-4,4]);
    yticks([-1:0.2:1.5]);
    annotation('textbox',[0.2, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
        'color',[0 0 0],'LineStyle','none','FontSize',18);
    annotation('textbox',[0.2, 0.2, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
        'color',[0 0 0],'LineStyle','none','FontSize',8)
    %     set(gcf,'Position',[100, 100, 900, 300]);
    
    str = [dir,...
        '1d_theta=',num2str(theta),...
        ];
    set(gca,'ydir','reverse')
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
    close all;
end


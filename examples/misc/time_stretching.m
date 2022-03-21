clear all;
close all;
dir = 'results/stretching_smooth_square/';
[~,~,~] = mkdir(dir);
ii = 1;

% x0 = [0.6343    0.2234   -1.5047    0.8396    0.5709   -0.7804    1.1095    2.3095    0.7147];
% x0 = [-1.1615    7.2705   -0.2378    0.1256    0.1357   0.2670    0.1610    2.9370    6.2626];
% x2 = [-7.0352   -2.2122    1.0642   -0.1090   -0.2427   0.3339   -0.1892   -3.0126   -1.1743];

% for gaussian
% x0  = [0.0461713906311539 0.0971317812358475 0.823457828327293 0.694828622975817];
% x1  = [-0.0461713906311539 0.0971317812358475 0.823457828327293 0.694828622975817];
% x2  = [0.0461713906311539 -0.0971317812358475 0.823457828327293 0.694828622975817];

% for square
% x0  =  [0 1 0.2 -0.1];
% x1  =  [0 -1 0.4 -0.1];
% x2  =  [0 1 0.2 0.1];

% for smooth square 
x0 = [0.3804 0.5678 0.0759, 0];
x1 = [0.3804 2 * 0.5678 0.0759, 0];
x2 = [0.3804 0.5678 -0.0759, 0];

% x0 = [-7.0352   -2.2122    1.0642   -0.1090   -0.2427   0.3339   -0.1892   -3.0126   -1.1743];
% x1 = x0;


theta_chart = [];
for angle = 0 : 2: 178
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
    
    offset = 100;
    theta_chart = vertcat(theta_chart,[mod(fix(2*(-offset + angle)),360),angle]);
    
end
theta_chart = sortrows(theta_chart);


for angle = [0:2:178]
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
    
    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);
    
    offset = 100;
    theta = mod(fix(2*(-offset + angle)),360);
    theta_str = num2str(theta);
    
    tc = time(time > -0.5 | time < 0.6);
    dc = dataset(time > -0.5 | time < 0.6,:);
    [me, ind] = max(dc,[],2) ;
    eind = EnergyCrop(ind);
    
    tc_all = time(time > -1 & time < 1.5);
    dc_all = dataset(time > -1 & time < 1.5,:);
    [me_all, ind_all] = max(dc_all,[],2) ;
    eind_all = EnergyCrop(ind_all);
    
    
    
    wb = abs(eind);
    
    [fitresult1, gof1] = createFit_smoothsquare(tc, eind, wb, x0);
    [fitresult2, gof2] = createFit_smoothsquare(tc, eind, wb, x1);
    [fitresult3, gof3] = createFit_smoothsquare(tc, eind, wb, x2);
    [fitresult4, gof4] = createFit_smoothsquare(tc, eind, wb, rand(1)*x0);
    
    
    if gof1.rmse < gof2.rmse
        gof = gof1;
        fitresult = fitresult1;
    else
        gof = gof2;
        fitresult = fitresult2;
    end
    
    if gof3.rmse < gof.rmse
        gof = gof3;
        fitresult = fitresult3;
    end
    
    if gof4.rmse < gof.rmse
        gof = gof4;
        fitresult = fitresult4;
    end
    %     [fitresult, gof] = createFit2(tc, eind, wb, x0);
    x0 = coeffvalues(fitresult);
    fitdata = fitresult(tc_all);
    
%     [fitresult, gof] = fit(tc(:), eind(:), 'poly1');
%     fitdata = fitresult(tc_all);
    
    close all;
    imagesc(EnergyCrop,time,dataset);
    annotation('textbox',[0.2, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',18);
    annotation('textbox',[0.2, 0.2, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',8);
    annotation('textbox',[0.6, 0.2, 0.1, 0.1],'String',['\sigma = ',num2str(fitresult.sigma,'%.2f')],...
        'color',[1 1 1],'LineStyle','none','FontSize',20)
    
    colormap('bone');
    hold on;
    
    plot( eind_all,tc_all,'LineStyle','none',...
        'Marker',"o",...
        'MarkerFaceColor',[0.8,0,0],...
        'MarkerSize',5,...
        'Color',[0.8,0,0]);
    plot(fitdata,tc_all,'LineWidth',2,'Color','#1a1aff');
    colorbar;
    set(gca,'FontSize',14);
    xlabel('Energy shift [eV]','FontSize',22);
    ylabel('Time delay [ps]','FontSize',22);
    yticks([-1 : 0.2 : 1.5]);
    
    xticks([-4:2:4]);
    ylim([-1,1.5]);
    xlim([-5,5]);
    %     set(gca,'Position',[100, 100, 500, 200]);
    hold off;
    
    str = [dir,...
        'theta=',num2str(theta),...
        ];
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
    
    res_struct(ii).('theta') = theta;
    res_struct(ii).('sigma') = fitresult.sigma;
    res_struct(ii).('rmse') = gof.rmse;
    ii = ii + 1;
    
end

save(strcat(dir,'res.mat'),"res_struct");

clear gca gcf;
dir = 'results/stretching_smooth_square/';

neg_dir = 'results/NEG_stretching_smooth_square/';
neg_res_struct = load(strcat(neg_dir,'res_NEG.mat')).('res_struct');
neg_res_struct = table2struct(sortrows(struct2table(neg_res_struct),'theta'));

load(strcat(dir,'res.mat'));

neg_angles = sort([[0 : 4: 28],[160 : 4: 200],[344 : 4: 356]].');
anglex_neg = [neg_res_struct(:).('theta')].';
sigma_neg = abs([neg_res_struct(:).('sigma')]).';
rmse_neg  = [neg_res_struct(:).('rmse')].';

rmse_neg = rmse_neg((anglex_neg >= 0 & anglex_neg <= 28) | ...
    (anglex_neg >= 160 & anglex_neg <= 200) | ...
    (anglex_neg >= 344 & anglex_neg <= 356));
sigma_neg = sigma_neg((anglex_neg >= 0 & anglex_neg <= 28) | ...
    (anglex_neg >= 160 & anglex_neg <= 200) | ...
    (anglex_neg >= 344 & anglex_neg <= 356));
% anglex_neg = anglex_neg((anglex_neg >= 0 & anglex_neg <= 28) | ...
%     (anglex_neg >= 160 & anglex_neg <= 200) | ...
%     (anglex_neg >= 344 & anglex_neg <= 356) );

figure;
res_struct = table2struct(sortrows(struct2table(res_struct),'theta'));
anglex = [res_struct(:).('theta')];


sigma = abs([res_struct(:).('sigma')]);
rmse  = [res_struct(:).('rmse')];

sigma((anglex_neg >= 0 & anglex_neg <= 28) | ...
    (anglex_neg >= 160 & anglex_neg <= 200) | ...
    (anglex_neg >= 344 & anglex_neg <= 356)) = sigma_neg;


rmse((anglex_neg >= 0 & anglex_neg <= 28) | ...
    (anglex_neg >= 160 & anglex_neg <= 200) | ...
    (anglex_neg >= 344 & anglex_neg <= 356)) = rmse_neg;

sigma2 = sigma;
anglex2 = anglex;

% anglex = anglex(sigma < 0.5 & sigma > 0);
% rmse = rmse(sigma < 0.5 & sigma > 0);
% sigma = sigma(sigma < 0.5 & sigma > 0);

yyaxis right;
stem(anglex,rmse,'LineWidth',0.5);
ylabel('Error (RMS)','FontSize',18);
xlabel('\theta [deg]','FontSize',18);
ylim([0,1]);
yticks([0:0.2:1]);
yyaxis left;

p = plot(anglex,sigma,'LineWidth',3,'Marker','*');
k = p.Parent;
k.Parent.Position = [360 198 800 300];
set(gca,'FontSize',14);
set(gca,'XMinorTick','on');

set(gca,'TickLength', [0.02, 0.06]);
ylabel('\sigma','FontSize',18);
ax = gca;
ax.XAxis.MinorTickValuesMode = 'manual';
ax.XAxis.MinorTickValues = 0:5:360;
xticks(unique(sort([[0 : 30: 360],[0:45:360]])));
xlim([0 , 360]);
ylim([0 , 0.6]);
exportgraphics(gcf, strcat(dir,'analysis','.png'),'resolution' , 400);
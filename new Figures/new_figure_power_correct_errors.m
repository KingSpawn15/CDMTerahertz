clear all 
close all
%%
delete(gcp('nocreate'));
parpool(6);

%%
% miscellaneous params
e_w = linspace(-5,5,181);

% photodember_parameters
pump_power_nj = 10;
laser_spot_size_fwhm = 40e-6;
fitting_parameter_EPD = (1.26/6.34) * 1.2;

% rectification parameters

params_rectification.tau = 30e-15;
params_rectification.lambda = 800e-9;
params_rectification.d  = .5e-3;
params_rectification.sigma_z = 50e-6;
params_rectification.z = -1e-6;
E_max_rectification = (1.631e6) * 1.65;
delay_or_pd_ps = 0.1;

% pulse energy list 
pulse_energy_list = [0.1000    0.1700    0.3000    0.5200    1.0000    1.7300    3.0000    5.1900   10.0000 17.3 30];

% Rectification field calculation
%%
[TOR, ZOR, EOR_10] = electric_field_rectification(params_rectification, ...
        E_max_rectification, delay_or_pd_ps);

% factor_mat = [1.5413    1.8102    2.5612    2.9008    1.7859    1.3690    1.4191    1.3701    1.2639    1.0813    0.9709];
factor_mat = [1.5413    1.8102    2.5612    2.9008    1.7859    1.3690    1.4191    1.3701    1.1    1.0813   0.9709 * 0.9];

for ii = 1 : length(pulse_energy_list)
    EOR{ii} = EOR_10 * (pulse_energy_list(ii) / 10)^(1/2) * (factor_mat(ii)) ;
%     EOR{ii} = EOR_10 * (pulse_energy_list(ii) / 10) / ((pulse_energy_list(ii) / 10) + ) ;
end

% 
% photodember field calculation
for ii = 1 : length(pulse_energy_list)
    eels_photodember = setup_parameters_eels_photodember(pulse_energy_list(ii), laser_spot_size_fwhm);
    [TPD, ZPD, EPD_store{ii}] = electric_field_photodember(eels_photodember, fitting_parameter_EPD);

end

%%
% photodember field calculation
for ii = 1 : length(pulse_energy_list)
    [T, Z, EPD{ii}, EOR{ii}] = interpolate_field(TOR, ZOR, EOR{ii}, TPD, ZPD, EPD_store{ii});
end


%% Get measurement data
[com,deltat,energy,eels_measure,errs] = getmeasurement_power_2();
[com_2,deltat_2,energy_2,eels_measure_2,errs_2] = getmeasurement_power();
%%

% [e_w, t_w, psi_incoherent_pd] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, 0 * EOR, T, Z, 45, eels_photodember, e_w);

for ii = 1:length(pulse_energy_list)
    disp(ii)
%     [TPD, ZPD, EPD{ii},psi_sub_pd{ii}, psi_incoherent_pd{ii}, ...
%         eels, w, e_w, t_w, tt, zz] = electric_field_photodember(pulse_energy_list(ii), ...
%         pd_spot_fwhm, pd_z_max);
%     EPD{ii} = EPD{ii}*(1.26/6.34);

[e_w, t_w, psi_incoherent_comb{ii}] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD{ii}, EOR{ii}, T, Z, 90, eels_photodember, e_w);
[t0_vec,eels_comb{ii}] = utils_spectrum.calculate_spectrum_from_fields(EOR{ii} + EPD{ii}, T, Z * 1e-6);

end


%%
[com_sim,errs_sim] = get_errors_sim(psi_incoherent_comb, e_w, t_w);
%%
close all
figure;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(3,4,"TileSpacing","compact");
% ax1 = axes(ttt);
% ax1.Layout.Tile = 5;
ttt.Padding = "loose";
nexttile
imagesc(energy, deltat - 0.3, eels_measure{5});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile;
imagesc(energy, deltat - 0.3, eels_measure{8});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile;
imagesc(energy, deltat- 0.3, eels_measure{9});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
% nexttile
% imagesc(energy, deltat- 0.3, eels_measure{10});
%     ylim([-1,1.5]);
%     xlim([-5,5]);
%     set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
%     colormap jet
%     axis square

nexttile
imagesc(energy, deltat- 0.3, eels_measure{11});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square

% nexttile(5,[3,3]);
% plot(pulse_energy_list, ecomb_max,'LineStyle','none', ...
%     'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
% hold on;
% fplot(@(x)saturate(param,x), [min(xdata), max(xdata)+5],LineWidth=2);
% set_axis_properties(gca,FontSize + 2,FontName,1,1e6*[0:0.5:6],0:5:30,'','',FontSize,[0.3 0.3 0.3])
% xlim([-5,35])

nexttile
imagesc(e_w,t_w, psi_incoherent_comb{5});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3])
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{8});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{9});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3]);

nexttile
imagesc(e_w,t_w, psi_incoherent_comb{11});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
axis square

nexttile([1,4]);
% errors = 0.1 * ones(size(deltat));
colors = turbo(11);

ERROR_VAL = 0;

ind = 1;
for i = [2     4     6     8     9    10  11]
%     plot(eels_comb{i} + 5 * (ind-1), t0_vec, 'color', colors(i,:),'LineWidth',1);
    
% errorbar(com(i,2:3:end) + 5 * (ind-1), t_w(2:3:end),errs_sim(i,2:3:end)/2, 'horizontal','color', ...
%         colors(i,:),'LineWidth',.5,'LineStyle','none');
%     hold on;
    
%     [errorx, errory] = errorband(com(i,:) + 5 *  (ind-1), errs(i,:), deltat - 0.3);

    [errorx, errory] = errorband(com_sim(i,:) + 5 *  (ind-1), errs_sim(i,:), t_w');
    fill(errorx, errory,  colors(i,:),'FaceAlpha',0.3,'EdgeColor','none');

    hold on;
    errorbar(com(i,2:2:end) + 5 *  (ind-1), deltat(2:2:end) - 0.3, ...
       ERROR_VAL + 1*errs(i,2:2:end)/2, 'horizontal',  'color', colors(i,:),'LineStyle','None','LineWidth',.5);

set(gca, 'YDir', 'reverse');
    ind = ind + 1
end
hold off
set_axis_properties(gca,FontSize,FontName,.2,-1:0.5:1.5,0:5:30,'','',FontSize,[0 0 0])
% axis equal;
% daspect([1 1 1]);
ylim([-1, 1]);
xlim([-4,35]);
set(gca,'XAxisLocation','top')


set(gcf,'Position',[200,50,200 + 600,200 + 600]); %set paper size (does not affect display)

exportgraphics(gcf, 'new Figures/results/power_fig_4.png', 'Resolution',300);
%%
% error analysis figure
close all
figure;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(1,2,"TileSpacing","compact");
% ax1 = axes(ttt);
% ax1.Layout.Tile = 5;
ttt.Padding = "loose";
nexttile
imagesc(energy, deltat- 0.3, eels_measure{10});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
hold on
errorbar(com_2(10,2:1:end),deltat(2:1:end) - .3,errs(10,2:1:end)/2,'horizontal', ...
    'Color',[1 1 1] * .7, ...
    'LineWidth',1.5,'LineStyle','none');

nexttile
imagesc(energy, deltat- 0.3, eels_measure{10});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
hold on
errorbar(com(10,2:1:end),deltat(2:1:end) - .3,errs(10,2:1:end)/2,'horizontal', ...
    'Color',[1 1 1] * .7, ...
    'LineWidth',1.5,'LineStyle','none');
set(gcf,'Position',[200,200,200 + 400,200 + 200]);
exportgraphics(gcf, 'article_check/results/error_analysis.png', 'Resolution',300);

max_exp = []
for i = 1:11
    max_exp = [max_exp, max(abs(com_2(i,:)))]
end
%%



%%
eor_max = zeros(11,1);
epd_max = zeros(11,1);
ecomb_max = zeros(11,1);
for i = 1:11
    Ecomb = EOR{i} + EPD{i};

    ecomb_max(i) =  max(abs(Ecomb(:)));

    eor_max(i) = max(abs(EOR{i}(:)));
    epd_max(i) = max(abs(EPD{i}(:)));

end

power = @(param,xdata) param(1).*(xdata.^param(2));
saturate = @(param,xdata) param(1).*(xdata)./(xdata + abs(param(2)));

% Initial guess for the parameters [A, Gamma, x0]
param0 = [1e7, 11];

% Your xdata and ydata here
xdata = pulse_energy_list(1:end).';
ydata = ecomb_max(1:end);

% Use lsqcurvefit to fit the Lorentzian to your data
options = optimset('Display','off');
[param,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(saturate, param0, xdata, ydata, [], [], options);


N = length(xdata) - length(param);
stdErrors = full(sqrt(diag(inv(J'*J)*resnorm/N)));

close all
figure;
plot(pulse_energy_list, ecomb_max,'LineStyle','none', ...
    'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
hold on;
fplot(@(x)saturate(param,x), [min(xdata), max(xdata)],LineWidth=2);
fplot(@(x) param(1) .*x ./ param(2), [min(xdata), max(xdata)],'--',LineWidth=2);
xlim([-5, 35]);
ylim([0, param(1) * 1.5]);
%%


close all;
figure;
T = tiledlayout(1,1,"TileSpacing","compact");
ax3 = axes(T);
ax1 = axes(T);
ax2 = axes(T);



plot(ax2,pulse_energy_list, eor_max,'LineStyle','none', ...
    'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','k'); 
hold on;


plot(ax2,pulse_energy_list, epd_max * 10,'LineStyle','none', ...
    'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k'); 
ax2.Box = 'off';
ax2.Color = 'none';
ax2.XTick = [];
ax2.YTick = [];
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';


fplot(ax1,@(x)saturate(param,x/159), [min(xdata) * 159, max(xdata)*159],'k--',LineWidth=2);
ax1.Box = 'off';
ax1.Color = 'none';
ax1.XTick = [];
ax1.YTick = [];
ax1.XAxis.Visible = 'off';
ax1.YAxis.Visible = 'off';


plot(ax3, pulse_energy_list, max_exp, 'Color','#0000CD','LineWidth',2, ...
    'MarkerSize',8,'Marker','s','MarkerEdgeColor','#00008B','MarkerFaceColor','#6495ED'); 
ax3.Box = 'off';
ax3.Color = 'none';
ax3.XTick = [];
ax3.YTick = [];
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'off';

ax2.XAxisLocation = 'bottom';
ax2.XAxis.Visible = 'on';
ax2.XTick = [0:5:30];
ax2.YLim = [0,7e6];

ax1.XAxisLocation = 'top';
ax1.XAxis.Visible = 'on';
ax1.XTick = [0:5:30] * 200;
ax1.YAxisLocation = 'left';
ax1.YAxis.Visible = 'on';
ax1.YTick = [0:1:7]*1e6;
ax1.YLim = [0,7e6];


ax3.YAxisLocation = 'right';
ax3.YAxis.Visible = 'on';
ax3.YTick = [0:0.5:5];
ax3.YLim = [0,5];
set(ax3,'YColor','#0000CD');

ax1.FontSize = FontSize+4;
ax2.FontSize = FontSize+4;
ax3.FontSize = FontSize+4;


AxisObject = ax1;   % Asigns the axis to a variable
ExponentToBeUsedLater = AxisObject.YAxis.Exponent; % Stores the exponent for later, since it is going to be erased
AxisObject.YTickLabelMode = 'manual'; 

set(gcf,'Position',[200,50,200 + 400,200 + 650]); %set paper size (does not affect display)
exportgraphics(gcf, 'new Figures/results/saturatin_fit.png', 'Resolution',300);
%%

% [com,deltat, energy, eels_measure,errs] = getmeasurement_power_2()

%%
function [errorx, errory] = errorband(mean, error, deltat)
    errorx = [mean + error/2, flip(mean - error/2)];
    errory = [deltat, flip(deltat)];
end


function [com,deltat, energy, eels_measure,errs] = getmeasurement_power()
    load('saved_matrices\PulseEnergy.mat')
    deltat = Time;
    energy = EnergyCrop;

    
    x_c = zeros(11, size(Time,2));
    errs = ones(11, size(Time,2));

    for ii = 1:11
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        eels_measure{ii} = eels;
         
        mass_matrix = eels;
        % Initialize a vector to store the x coordinates of the center of mass of each row
%         errs = zeros(1, size(deltat,2));
        for i = 1:size(eels,1)
%             row_mass = sum(mass_matrix(i,:));
            for j = 1:size(eels,2)
                [~, ind_max] = max(mass_matrix(i,:));
                x_c(ii,i) = ind_max;
            end
            x_c(ii,i) = x_c(ii,i);

            errs(ii,i) = energy(find(mass_matrix(i,:) > 0.7 * max(mass_matrix(i,:)),1,'last')) - ...
            energy(find(mass_matrix(i,:) >  0.7 * max(mass_matrix(i,:)),1,'first'));

        end
    end

    com = energy(floor(x_c));
end


function [com,errs] = get_errors_sim(eels_sim, e_w, t_w)
        
    com = zeros(11, size(t_w,1));
    errs = ones(11, size(t_w,1));
    
    for ii = 1:11
        eels = eels_sim{ii};
        eels = eels ./ sqrt(sum(eels.^2, 2));
        [errs_i,means_i] = error_calculator(eels,e_w);

        com(ii,:) = means_i;
        errs(ii,:) = errs_i;
    
    end

end




function [com,deltat, energy, eels_measure,errs] = getmeasurement_power_2()
    load('saved_matrices\PulseEnergy.mat')
    
    deltat = Time;
    energy = EnergyCrop;
    
    com = zeros(11, size(Time,2));
    errs = ones(11, size(Time,2));

    for ii = 1:11
        
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        [errs_i,means_i] = error_calculator(eels,energy);

        eels_measure{ii} = eels;
        com(ii,:) = means_i;
        errs(ii,:) = errs_i;
    
    end
% 
%     [errs,means] = error_calculator(psi_in,e_w)
%     
%     x_c = zeros(11, size(Time,2));
%     errs = ones(11, size(Time,2));
% 
%     for ii = 1:11
%         eels = squeeze(cell2mat(DataSetCropAll(ii)))';
%         eels = eels ./ sqrt(sum(eels.^2, 2));
%         eels_measure{ii} = eels;
%          
%         
% 
%         end
%     end
% 
%     com = energy(floor(x_c));
end

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end

function [errs,means] = error_calculator(psi_in,e_w)
%ERROR_CALCULATOR Summary of this function goes here
%   Detailed explanation goes here
errs = []
means = []

for ind = 1 : size(psi_in,1)
    psi_in(psi_in<0) = 0;
    psi = psi_in(ind, :);
    errs = [errs,std(e_w, psi / sum(psi))];
    means = [means,dot(e_w , psi /sum(psi))];
end

end


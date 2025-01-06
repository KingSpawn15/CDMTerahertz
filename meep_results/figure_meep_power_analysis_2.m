clearvars
close all


params.e_w = linspace(-5, 5, 271);
params.exp_theory_time_shift = 0.1;

common_fac = 1.46 * 0.9281;
optimal_parameters.weight_pd = common_fac * (-0.09 * 1e18);
optimal_parameters.weight_or = common_fac * (1 * 1e9);

[~, ~, EPD_xz] = get_fields_photodember_meep();

e_w = params.e_w;

% Parameter ranges
shift_values = -32;
adjust_or_values = 1;1.2;
adjust_pd_values = 0.7;0.8;
spotsize_values = 117;


setdir = 'meep_results/results/adjustments/';


params.spot_size_fwhm_um_or = spotsize_values;
[tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(spotsize_values);
[T, Z] = ndgrid(tc, xc);


close all;
% Apply adjustments
EPD_s = shiftWithZeros2D(EPD_xz * optimal_parameters.weight_pd * adjust_pd_values, 0, shift_values);
EOR_s = (EOR_zz) * optimal_parameters.weight_or * adjust_or_values;

%%
pulse_energy_list = [0.10    0.17    0.30    0.52    1.00    1.73    3.00    5.19   10.00 17.3 30.0];
max_epd = zeros(1,11);
max_eor = zeros(1,11);
max_etot = zeros(1,11);

ps_or = 4.7160;
ps_pd = 3;
% A = (10 + ps) / 10;
weight = @(pow, ps) ((10 + ps) / 10) * pow / (pow + ps);

psi_incoherent_comb= cell(1, 11); 
EPD = cell(1, 11); 
EOR = cell(1, 11); 

for ii = 1:length(pulse_energy_list)

    EPD{ii} = EPD_s * weight(pulse_energy_list(ii), ps_pd);
    EOR{ii} = EOR_s * weight(pulse_energy_list(ii), ps_or);

    [t_w_0 , psi_incoherent_comb{ii}] = calculate_incoherent_spectrum_from_fields(EOR{ii} + EPD{ii}, T, Z, e_w);
    max_epd(ii) = max(abs(EPD{ii}(:)));
    max_eor(ii) = max(abs(EOR{ii}(:)));
    max_etot(ii) = max(abs(EPD{ii}(:) + EOR{ii}(:)));
end


%% Fluence conversion factor for our beam.
energy_to_fluence_factor = 159.2;

%% Get measurement data

[com,deltat,energy,eels_measure,errs] = getmeasurement_power();
[com_m,deltat_m,energy_m,eels_measure_m,errs_m] = getmeasurement_power_max();

max_exp = zeros(11,1);

for ind = 1:11
    max_exp(ind) = max(com_m(ind,2:1:end));
end

saturate = @(param,xdata) param(1).*(xdata)./(xdata + abs(param(2)));

% Initial guess for the parameters [A, Gamma, x0]
param0 = [1, 5];


xdata = pulse_energy_list(1:end).';
ydata = max_exp(1:end);

options = optimset('Display','off');
[param,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(saturate, param0, xdata, ydata, [], [], options);

N = length(xdata) - length(param);
stdErrors = full(sqrt(diag(inv(J'*J)) * resnorm/N));


disp(['saturation fluence = ', num2str(param(2)  * energy_to_fluence_factor), '$\mu J / cm^2$'])
disp(['saturation fluence (err)= ', num2str(stdErrors(2)  * energy_to_fluence_factor), '$\mu J / cm^2$'])
%% sim errors if wish
[com_sim,errs_sim] = get_errors_sim(psi_incoherent_comb, e_w, t_w_0);

%% Figure 4(A,B,C)
t_w = t_w_0 - params.exp_theory_time_shift;
close all
figure;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(3,4,"TileSpacing","compact");
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
nexttile
imagesc(energy, deltat- 0.3, eels_measure{11});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
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
colors = turbo(11);

ERROR_VAL = 0;

ind = 1;
for i = [2     4     6     8     9    10  11]

    [errorx, errory] = errorband(com_sim(i,:) + 5 *  (ind-1), errs_sim(i,:), t_w' - 0.3);
    fill(errorx, errory,  colors(i,:),'FaceAlpha',0.3,'EdgeColor','none');

    hold on;
    errorbar(com(i,2:2:end) + 5 *  (ind-1), deltat(2:2:end) - 0.3, ...
       ERROR_VAL + 1*errs(i,2:2:end)/2, 'horizontal',  'color', colors(i,:),'LineStyle','None','LineWidth',.5);

    set(gca, 'YDir', 'reverse');
    ind = ind + 1;
end

hold off
set_axis_properties(gca,FontSize,FontName,.2,-1:0.5:1.5,0:5:30,'','',FontSize,[0 0 0])

ylim([-1, 1]);
xlim([-4,35]);
set(gca,'XAxisLocation','top')


set(gcf,'Position',[200,50,200 + 600,200 + 600]); %set paper size (does not affect display)
exportgraphics(gcf, 'meep_results/results/diffusion_model_power_fig_4.png', 'Resolution',300);
%% %% Figure 4(D)

close all;
figure;
T = tiledlayout(1,1,"TileSpacing","compact");
ax3 = axes(T);
hold(ax3,'on');
ax1 = axes(T);
ax2 = axes(T);


plot(ax2,pulse_energy_list, max_eor,'LineStyle','none', ...
    'Marker','d','MarkerEdgeColor','#32CD32','MarkerFaceColor','#FFFFFF','MarkerSize',12, 'LineWidth', 2); 
hold on;
plot(ax2,pulse_energy_list, max_epd * 50,'LineStyle','none', ...
    'Marker','o','MarkerEdgeColor','#707070','MarkerFaceColor','#707070','MarkerSize',10); 
ax2.Box = 'off';
ax2.Color = 'none';
ax2.XTick = [];
ax2.YTick = [];
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';


fplot(ax1,@(x)saturate(param,x), [min(xdata) * energy_to_fluence_factor, max(xdata) * energy_to_fluence_factor],'r--','LineStyle','none');
ax1.Box = 'off';
ax1.Color = 'none';
ax1.XTick = [];
ax1.YTick = [];
ax1.XAxis.Visible = 'off';
ax1.YAxis.Visible = 'off';

xx = 0:.5:35;
yy = saturate(param,xx);
plot(ax3, pulse_energy_list, max_exp, 'Color','#0000CD','LineStyle','None', ...
    'MarkerSize',14,'Marker','s','MarkerEdgeColor','#0000CD','MarkerFaceColor','#0000CD'); 
plot(ax3, xx, yy, 'Color','#0000CD','LineWidth',2,'LineStyle','-'); 
ax3.Box = 'off';
ax3.Color = 'none';
ax3.XTick = [];
ax3.YTick = [];
ax3.XLim = [0,35];
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'off';

ax2.XAxisLocation = 'bottom';
ax2.XAxis.Visible = 'on';
ax2.XTick = 0:5:30;
ax2.YLim = [0,25e6];
ax2.YTick = (0:5:25)*1e6;
ax2.XLim = [0,35];
ax2.YAxis.Visible = 'on';

ax1.XAxisLocation = 'top';
ax1.XAxis.Visible = 'on';
ax1.XTick = (0:5:30) * 200;
ax1.YAxisLocation = 'left';
ax1.YAxis.Visible = 'off';
ax1.YTick = (0:1:7)*1e6;
ax1.YLim = [0,7e6];


ax3.YAxisLocation = 'right';
ax3.YAxis.Visible = 'on';
ax3.YTick = 0:0.5:5;
ax3.YLim = [0,5];
set(ax3,'YColor','#0000CD');

grid(ax3, 'on');
ax3.GridColor = [0.8 0.8 0.8]; % Light gray grid lines
grid(ax2, 'on');
ax2.GridColor = [0.8 0.8 0.8]; % Light gray grid lines

ax1.FontSize = FontSize+4;
ax2.FontSize = FontSize+4;
ax3.FontSize = FontSize+4;


AxisObject = ax2;   % Asigns the axis to a variable
ExponentToBeUsedLater = AxisObject.YAxis.Exponent; % Stores the exponent for later, since it is going to be erased
AxisObject.YTickLabelMode = 'manual'; 

set(gcf,'Position',[200,50,200 + 400,200 + 650]); %set paper size (does not affect display)
exportgraphics(gcf, 'meep_results/results/diffusion_model_saturatin_fit.png', 'Resolution',300);



%% error analysis figure
% 
% for ii = 1:11
%     [com_sim_max, errs_sim_max] = getsimulation_power_max(psi_incoherent_comb{ii}, e_w, t_w);
%     cs{ii} = com_sim_max ;
%     es{ii} = errs_sim_max;
% end
% 
% com_max_s = max(cell2mat(cs')');
% com_meas_s = max(com_m');
% close all
% 
% [com_sim_max, errs_sim_max] = getsimulation_power_max(psi_incoherent_comb{10}, e_w, t_w)
% 
% figure;
% set(groot,'defaultAxesXTickLabelRotationMode','manual')
% FontName = 'ariel';
% FontSize = 14;
% ttt = tiledlayout(1,2,"TileSpacing","compact");
% % ax1 = axes(ttt);
% % ax1.Layout.Tile = 5;
% ttt.Padding = "loose";
% nexttile
% imagesc(energy, deltat- 0.3, eels_measure{10});
%     ylim([-1,1.5]);
%     xlim([-5,5]);
%     set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3])
%     colormap jet
%     axis square
% hold on
% errorbar(com_m(10,2:1:end),deltat(2:1:end) - .3,errs_m(10,2:1:end)/10,'horizontal', ...
%     'Color',[1 1 1] * .7, ...
%     'LineWidth',1.5,'LineStyle','none');
% set(gcf,'Position',[200,200,200 + 400,200 + 200]);
% 
% nexttile
% imagesc(e_w, t_w, psi_incoherent_comb{10});
%     ylim([-1,1.5]);
%     xlim([-5,5]);
%     set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3])
%     colormap jet
%     axis square
% hold on
% errorbar(com_sim_max(2:1:end),t_w(2:1:end) ,errs_sim_max(2:1:end)/10,'horizontal', ...
%     'Color',[1 1 1] * .7, ...
%     'LineWidth',1.5,'LineStyle','none');

% exportgraphics(gcf, 'meep_results/results/error_analysis.png', 'Resolution',300);

% function [tc, xc, field_pd] = get_fields_photodember_meep_intensities()
%     
%     pulse_energy_list = [0.10    0.17    0.30    0.52    1.00    1.73    3.00    5.19   10.00 17.30 30.00];
%     formatSpec = '%.2f';
% %     dir = 'meep_results/saved_matrices_meep/photodember/test/'
%     dir = 'meep_results\saved_matrices_meep\photodember\combined\';
%     for pe = 1:length(pulse_energy_list)
% 
%         filename = strcat(dir,'field_ez_pd_intensity_',num2str(pulse_energy_list(pe),formatSpec),'t0_0.6fwhm_t_50.mat');
%         load(filename);
%         field_pd{pe} = e_pd.';
%     end
% 
%     xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
%     tc = tstep : tstep : tstep * size(e_pd,1);
% 
% 
% end

% function plot_tile(x, y, z)
%     FontSize = 18;
%     FontName = 'ariel';
%     nexttile
%     imagesc(x, y, z);
%     ylim([-1,1.5]);
%     xlim([-5,5]);
%     colormap jet
%     axis square
%     set(groot,'defaultAxesXTickLabelRotationMode','manual');
%     set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3]);
%     
% end

function set_axis_properties(ax,FontSize,FontName,~,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end

function [com,deltat, energy, eels_measure,errs] = getmeasurement_power()
    load('saved_matrices\PulseEnergy.mat', 'DataSetCropAll', 'Time', 'EnergyCrop');
    
    deltat = Time;
    energy = EnergyCrop;
    
    com = zeros(11, size(Time,2));
    errs = ones(11, size(Time,2));
    
    eels_measure = cell(1, 11); 

    for ii = 1:11
        
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        [errs_i,means_i] = error_calculator(eels,energy);

        eels_measure{ii} = eels;
        com(ii,:) = means_i;
        errs(ii,:) = errs_i;
    
    end

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


function [errorx, errory] = errorband(mean, error, deltat)
    errorx = [mean + error/2, flip(mean - error/2)];
    errory = [deltat, flip(deltat)];
end


function [errs,means] = error_calculator(psi_in,e_w)
    %ERROR_CALCULATOR Summary of this function goes here
    %   Detailed explanation goes here
%     errs = [];
%     means = [];
    
    % Preallocate arrays
    num_rows = size(psi_in, 1);
    errs = zeros(1, num_rows);
    means = zeros(1, num_rows);

    for ind = 1 : size(psi_in,1)
        psi_in(psi_in<0) = 0;
        psi = psi_in(ind, :);

        errs(ind) = std(e_w, psi / sum(psi));  % Store standard deviation
        means(ind) = dot(e_w, psi / sum(psi)); % Store weighted mean

    end

end

function [com,deltat, energy, eels_measure,errs] = getmeasurement_power_max()
    load('saved_matrices\PulseEnergy.mat', 'DataSetCropAll', 'Time', 'EnergyCrop')
    deltat = Time;
    energy = EnergyCrop;

    
    x_c = zeros(11, size(Time,2));
    errs = ones(11, size(Time,2));
    
    eels_measure = cell(1, 11); 

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

function [com, errs] = getsimulation_power_max(eels, e_w, t_w)

    energy = e_w;

    
    x_c = zeros(1, size(t_w,2));
    errs = ones(1, size(t_w,2));

        mass_matrix = eels;

        for i = 1:size(eels,1)
            for j = 1:size(eels,2)
                [~, ind_max] = max(mass_matrix(i,:));
                x_c(1,i) = ind_max;
            end
            x_c(1,i) = x_c(1,i);

            errs(1,i) = energy(find(mass_matrix(i,:) > 0.7 * max(mass_matrix(i,:)),1,'last')) - ...
            energy(find(mass_matrix(i,:) >  0.7 * max(mass_matrix(i,:)),1,'first'));

        end

    com = energy(floor(x_c));
end

function [tc, xc, field_pd] = get_fields_photodember_meep()

    load('meep_results\saved_matrices_meep\photodember\eps_12\field_ez_pd_intensity_x_10t0_0.5fwhm_t_50sigma_t_40.mat')
%     load('meep_results\saved_matrices_meep\photodember\test\field_ez_pd_intensity_10.00t0_0.5fwhm_t_50.mat')
%     load('meep_results\saved_matrices_meep\photodember\varying_intensities\field_ez_pd_intensity_5.19t0_0.5.mat')
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

end

function shiftedMatrix = shiftWithZeros2D(matrix, shiftRow, shiftCol)
    % Get the size of the matrix
    [rows, cols] = size(matrix);

    % Initialize the output matrix with zeros
    shiftedMatrix = zeros(rows, cols);

    % Calculate the ranges for rows and columns in the shifted matrix
    rowStart = max(1, 1 + shiftRow);
    rowEnd = min(rows, rows + shiftRow);
    colStart = max(1, 1 + shiftCol);
    colEnd = min(cols, cols + shiftCol);

    % Calculate the ranges for rows and columns in the original matrix
    originalRowStart = max(1, 1 - shiftRow);
    originalRowEnd = min(rows, rows - shiftRow);
    originalColStart = max(1, 1 - shiftCol);
    originalColEnd = min(cols, cols - shiftCol);

    % Copy the values from the original matrix to the shifted matrix
    shiftedMatrix(rowStart:rowEnd, colStart:colEnd) = ...
        matrix(originalRowStart:originalRowEnd, originalColStart:originalColEnd);
end

clear all
close all

params.e_w = linspace(-5,5,181);
params.exp_theory_time_shift = 0.6;
params.spot_size_fwhm_um_or = 70;

common_fac = 0.9;
optimal_parameters.weight_pd = 1.2 * 1 * common_fac ;
optimal_parameters.weight_or = 1e4 * .8 * 4 * 2. * common_fac ;

[tc, xc, EPD_xz] = get_fields_photodember_meep_intensities();
% [~, ~, EPD_xz] = get_fields_photodember_meep();
[~, ~, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(params.spot_size_fwhm_um_or);
[T, Z] = ndgrid(tc, xc);
e_w = params.e_w;
%%
pulse_energy_list = [0.10    0.17    0.30    0.52    1.00    1.73    3.00    5.19   10.00 17.3 30.0];
factor_mat = [1.5413    1.8102    2.5612    2.9008    1.7859    1.3690    1.4191    1.3701    1.0    1.1 * 10/17.3   1.12 * 10/30];

max_epd = []
max_eor = []
max_etot = []

for ii = 1:length(pulse_energy_list)
    intensity_weight_or = (pulse_energy_list(ii) / 10) * factor_mat(ii);
    EOR = EOR_zz * optimal_parameters.weight_or * intensity_weight_or;
    EPD = movmean(EPD_xz{ii},10,1) * optimal_parameters.weight_pd;
%     EPD(abs(Z')>40) = 0;
    [t_w_0 , psi_incoherent_comb{ii}] = calculate_incoherent_spectrum_from_fields(EOR + EPD, T, Z, e_w);
    max_epd = [max_epd , max(abs(EPD(:)))]
    max_eor = [max_eor , max(abs(EOR(:)))]
    max_etot = [max_etot , max(abs(EPD(:) + EOR(:)))]
end





% [t_w_0 , psi_incoherent_comb_0] = calculate_incoherent_spectrum_from_fields(-EOR + EPD, T, Z, e_w);
% [~ , psi_incoherent_comb_45] = calculate_incoherent_spectrum_from_fields(EPD, T, Z, e_w);
% [~ , psi_incoherent_comb_90] = calculate_incoherent_spectrum_from_fields(EOR + EPD, T, Z, e_w);
% 
% [~ , psi_incoherent_or_0] = calculate_incoherent_spectrum_from_fields(-EOR, T, Z, e_w);
% [~ , psi_incoherent_or_45] = calculate_incoherent_spectrum_from_fields(-EOR * 0, T, Z, e_w);
% [~ , psi_incoherent_or_90] = calculate_incoherent_spectrum_from_fields(EOR, T, Z, e_w);
% 
% [~ , psi_incoherent_pd] = calculate_incoherent_spectrum_from_fields(EPD, T, Z, e_w);

%% Get measurement data

[com,deltat,energy,eels_measure,errs] = getmeasurement_power();
%%
[com_sim,errs_sim] = get_errors_sim(psi_incoherent_comb, e_w, t_w_0 );

%%

t_w = t_w_0 - params.exp_theory_time_shift;
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
imagesc(energy, deltat - 0.3, eels_measure{8});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile;
imagesc(energy, deltat - 0.3, eels_measure{9});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile;
imagesc(energy, deltat- 0.3, eels_measure{10});
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
imagesc(e_w,t_w, psi_incoherent_comb{8});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3])
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{9});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{10});
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

    [errorx, errory] = errorband(com_sim(i,:) + 5 *  (ind-1), errs_sim(i,:), t_w');
    fill(errorx, errory,  colors(i,:),'FaceAlpha',0.3,'EdgeColor','none');

    hold on;
    errorbar(com(i,2:2:end) + 5 *  (ind-1), deltat(2:2:end) - 0.3, ...
       ERROR_VAL + 1*errs(i,2:2:end)/2, 'horizontal',  'color', colors(i,:),'LineStyle','None','LineWidth',.5);

    set(gca, 'YDir', 'reverse');
    ind = ind + 1;
end

hold off
set_axis_properties(gca,FontSize,FontName,.2,-1:0.5:1.5,0:5:30,'','',FontSize,[0 0 0])
% axis equal;
% daspect([1 1 1]);
ylim([-1, 1]);
xlim([-4,35]);
set(gca,'XAxisLocation','top')


set(gcf,'Position',[200,50,200 + 600,200 + 600]); %set paper size (does not affect display)

exportgraphics(gcf, 'meep_results/results/power_fig_4.png', 'Resolution',300);
%%
% close all
% imagesc(e_w, t_w_0, psi_incoherent_comb{11});
%     ylim([-1,1.5]);
%     xlim([-5,5]);
%     set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3])
%     colormap jet
%     axis square
% hold on
% errorbar(com_sim(11,2:1:end),t_w_0(2:1:end),errs_sim(10,2:1:end)/2,'horizontal', ...
%     'Color',[1 1 1] * .7, ...
%     'LineWidth',1.5,'LineStyle','none');
% set(gcf,'Position',[200,200,200 + 200,200 + 200]);


%%
% function [tc, xc, field_pd] = get_fields_photodember_meep()
% 
%     load('meep_results/saved_matrices_meep/photodember/spot_size_30_shift04/field_ez_pd_intensity_10t0_0.5.mat')
%     xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
%     tc = tstep : tstep : tstep * size(e_pd,1);
%     field_pd = e_pd.';
% 
% end

function [tc, xc, field_pd] = get_fields_photodember_meep_intensities()
    
    pulse_energy_list = [0.10    0.17    0.30    0.52    1.00    1.73    3.00    5.19   10.00 17.30 30.00];
    formatSpec = '%.2f'
%     dir = 'meep_results/saved_matrices_meep/photodember/test/'
    dir = 'meep_results\saved_matrices_meep\photodember\combined\';
    for pe = 1:length(pulse_energy_list)

        filename = strcat(dir,'field_ez_pd_intensity_',num2str(pulse_energy_list(pe),formatSpec),'t0_0.6fwhm_t_50.mat');
        load(filename);
        field_pd{pe} = e_pd.'
    end

    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);


end

function plot_tile(x, y, z)
    FontSize = 18;
    FontName = 'ariel';
    nexttile
    imagesc(x, y, z);
    ylim([-1,1.5]);
    xlim([-5,5]);
    colormap jet
    axis square
    set(groot,'defaultAxesXTickLabelRotationMode','manual');
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3]);
    
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

function [com,deltat, energy, eels_measure,errs] = getmeasurement_power()
    load('saved_matrices\PulseEnergy.mat');
    
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

%%
function [errorx, errory] = errorband(mean, error, deltat)
    errorx = [mean + error/2, flip(mean - error/2)];
    errory = [deltat, flip(deltat)];
end

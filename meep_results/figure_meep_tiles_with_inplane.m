clear all
close all

params.e_w = linspace(-5,5,181);
params.exp_theory_time_shift = 0.7;
params.spot_size_fwhm_um_or = 70;

common_fac = .9 * 3e2;
optimal_parameters.weight_pd = 1.2 * 1 * common_fac ;
optimal_parameters.weight_or = 1e4 * .8 * 4 * 2. * common_fac ;

[~, ~, EPD_xz] = get_fields_photodember_meep();
[tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(params.spot_size_fwhm_um_or);
[T, Z] = ndgrid(tc, xc);
e_w = params.e_w;
%%

EOR = (EOR_zz ) * optimal_parameters.weight_or;
EOR_ip = - (EOR_xz ) * optimal_parameters.weight_or * (1/sqrt(2));
EPD = EPD_xz * optimal_parameters.weight_pd;

[t_w_0 , psi_incoherent_comb_0] = calculate_incoherent_spectrum_from_fields(-EOR + EPD + EOR_ip, T, Z, e_w);
[~ , psi_incoherent_comb_45] = calculate_incoherent_spectrum_from_fields(EPD + EOR_ip, T, Z, e_w);
[~ , psi_incoherent_comb_90] = calculate_incoherent_spectrum_from_fields(EOR + EPD + EOR_ip, T, Z, e_w);

[~ , psi_incoherent_or_0] = calculate_incoherent_spectrum_from_fields(-EOR + EOR_ip, T, Z, e_w);
[~ , psi_incoherent_or_45] = calculate_incoherent_spectrum_from_fields(-EOR * 0 + EOR_ip, T, Z, e_w);
[~ , psi_incoherent_or_90] = calculate_incoherent_spectrum_from_fields(EOR + EOR_ip, T, Z, e_w);

[~ , psi_incoherent_pd] = calculate_incoherent_spectrum_from_fields(EPD, T, Z, e_w);

%% Get measurement data

[~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
[~, ~, ~, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));
[~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));


%%
    
t_w = t_w_0 - params.exp_theory_time_shift;
setdir = 'meep_results/results/';
close all
close all
figure;
image_name = strcat('tiles_meep_inplane_2');
FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(3,4,"TileSpacing","compact");
ttt.Padding = "loose";

plot_tile(energy, deltat, eels_measure_0);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3]);
plot_tile(e_w,t_w, psi_incoherent_comb_0);
plot_tile(e_w,t_w, psi_incoherent_pd);
plot_tile(e_w,t_w, psi_incoherent_or_0);

plot_tile(energy, deltat, eels_measure_45);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3]);
plot_tile(e_w,t_w, psi_incoherent_comb_45);
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
plot_tile(e_w,t_w, psi_incoherent_pd);
plot_tile(e_w,t_w, psi_incoherent_or_45);

plot_tile(energy, deltat, eels_measure_90);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
plot_tile(e_w,t_w, psi_incoherent_comb_90);
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
plot_tile(e_w,t_w, psi_incoherent_pd);
plot_tile(e_w,t_w, psi_incoherent_or_90);

set(gcf,'Position',[200,200,200 + 3 * 300,200 +  400]); %set paper size (does not affect display)


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

%%
function [tc, xc, field_pd] = get_fields_photodember_meep()

    load('meep_results\saved_matrices_meep\photodember\combined\field_ez_pd_intensity_10.00t0_0.6fwhm_t_50.mat')
%     load('meep_results\saved_matrices_meep\photodember\test\field_ez_pd_intensity_10.00t0_0.5fwhm_t_50.mat')
%     load('meep_results\saved_matrices_meep\photodember\varying_intensities\field_ez_pd_intensity_5.19t0_0.5.mat')
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

end

% function [tc, xc, field_pd] = get_fields_photodember_meep_intensities()
% 
%     load('meep_results/saved_matrices_meep/photodember/spot_size_30_shift04/field_ez_pd_intensity_10t0_0.5.mat')
%     xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
%     tc = tstep : tstep : tstep * size(e_pd,1);
%     field_pd = e_pd.';
% 
% end

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


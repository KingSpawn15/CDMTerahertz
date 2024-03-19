clear all
close all

params.e_w = linspace(-5,5,181);
params.exp_theory_time_shift = 0.6;
params.spot_size_fwhm_um_or = 80;

optimal_parameters.weight_pd = 1.2;
optimal_parameters.weight_or = 1e3 * .8;

[~, ~, EPD_xz] = get_fields_photodember_meep();
[tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(params.spot_size_fwhm_um_or);
[T, Z] = ndgrid(tc, xc);
e_w = params.e_w;


EOR = EOR_zz * optimal_parameters.weight_or;
EPD = EPD_xz * optimal_parameters.weight_pd;



%% Get measurement data

angle = sort([0, 30, 50, 70, 90, 110, 130, 150, 180]);

% [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
for ii = 1:length(angle)
[~, deltat, energy, eels_measure{ii}] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(angle(ii)));
end


%%

for ii = 1:length(angle)
    disp(angle(ii))
    [t_w_0 , psi_incoherent_comb{ii}] = calculate_incoherent_spectrum_from_fields(-EOR * cos(2 * angle(ii) * pi /180) + EPD, T, Z, e_w);
end


%% Get measurement data

[~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
[~, ~, ~, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));
[~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));


%%
    
t_w = t_w_0 - params.exp_theory_time_shift;
setdir = 'meep_results/results/';
figure;
image_name = 'tiles_multiple_meep';
FontName = 'ariel';
FontSize = 15;
ttt = tiledlayout(2,length(angle),"TileSpacing","compact");
ttt.Padding = "compact";

for ii = 1:length(angle)
plot_tile(energy, deltat, eels_measure{ii});

set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[],'','',FontSize,[0 0 0, 0]);
end

end

for ii = 1:length(angle)
plot_tile(e_w,t_w, psi_incoherent_comb{ii});
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[-4:2:4],'','',FontSize,[0 0 0, 0]);
end
end


set(gcf,'Position',[200,200,200 + 120 * length(angle), 200 +  160]);


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',500);

%%
function [tc, xc, field_pd] = get_fields_photodember_meep()

    load('meep_results/saved_matrices_meep/photodember/spot_size_30_shift04/field_ez_pd_intensity_10t0_0.5.mat')
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

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


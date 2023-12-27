clear all
close all

fields = optimal_parameters();
EOR = fields.EOR;
EPD = fields.EPD;
e_w = fields.e_w ;
T = fields.T ;
Z = fields.Z ;
eels_photodember = fields.eels_obj;
%% Get measurement data

[~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
[~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));


%%
[e_w, t_w, psi_incoherent_pd] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, 0 * EOR, T, Z, 45, eels_photodember, e_w);
[~, ~, psi_incoherent_or_0] = utils_spectrum.generate_incoherent_spectrum_for_angle(0 * EPD, EOR, T, Z, 0, eels_photodember, e_w);
[~, ~, psi_incoherent_or_90] = utils_spectrum.generate_incoherent_spectrum_for_angle(0 * EPD, EOR, T, Z, 90, eels_photodember, e_w);
[~, ~, psi_incoherent_comb_0] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, EOR, T, Z, 0, eels_photodember, e_w);
[~, ~, psi_incoherent_comb_90] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, EOR, T, Z, 90, eels_photodember, e_w);

%%
setdir = 'new Figures/results/';
close all
close all
figure;
image_name = 'tiles';
FontName = 'ariel';
FontSize = 16;
ttt = tiledlayout(2,4,"TileSpacing","compact");
ttt.Padding = "loose";

plot_tile(energy, deltat, eels_measure_0);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3]);
plot_tile(e_w,t_w, psi_incoherent_comb_0);
plot_tile(e_w,t_w, psi_incoherent_pd);
plot_tile(e_w,t_w, psi_incoherent_or_0);
plot_tile(energy, deltat, eels_measure_90);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
plot_tile(e_w,t_w, psi_incoherent_comb_90);
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
plot_tile(e_w,t_w, psi_incoherent_pd);
plot_tile(e_w,t_w, psi_incoherent_or_90);

set(gcf,'Position',[200,200,200 + 2 * 300,200 +  300]); %set paper size (does not affect display)


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

%%
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


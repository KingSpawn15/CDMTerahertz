clearvars -except TPD_non_interpolated ZPD_non_interpolated EPD_non_interpolated eels_photodember
close all

%  photodember_parameters
% delete(gcp('nocreate'));
% parpool(28)
% pump_power_nj = 10;
% laser_spot_size_fwhm = 40e-6;
% fitting_parameter_EPD = (1.26/6.34) * 1.2;
% eels_photodember = setup_parameters_eels_photodember(pump_power_nj, laser_spot_size_fwhm);
% [TPD_non_interpolated, ZPD_non_interpolated, EPD_non_interpolated] = electric_field_photodember(eels_photodember, fitting_parameter_EPD);
%%
spot_size = 80;
shift = 0.58 ;
weight = 1 * (30/50)^(5/2);
fields = optimal_parameters_dual(weight, spot_size, shift, TPD_non_interpolated, ZPD_non_interpolated, ...
    EPD_non_interpolated, eels_photodember);
EOR_0 = -fields.EOR_90 ;
EOR_90 = fields.EOR_90 ;
EPD = fields.EPD * 1.2;
e_w = fields.e_w ;
T = fields.T ;
Z = fields.Z ;
eels_photodember = fields.eels_obj;
%%

%% Get measurement data

[~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
[~, ~, ~, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));
[~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));

%%
EOR_store_0 = EOR_0;
EOR_store_90 = EOR_90;
%%
EOR_0 = EOR_store_0 * .5;
EOR_90 = EOR_store_90 * .5;

[e_w, t_w, psi_incoherent_pd] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, 0 * EOR_0, T, Z, 45, eels_photodember, e_w);
[~, ~, psi_incoherent_or_0] = utils_spectrum.generate_incoherent_spectrum_for_angle(0 * EPD, EOR_0, T, Z, 0, eels_photodember, e_w);
[~, ~, psi_incoherent_or_45] = utils_spectrum.generate_incoherent_spectrum_for_angle(0 * EPD, EOR_0, T, Z, 45, eels_photodember, e_w);
[~, ~, psi_incoherent_or_90] = utils_spectrum.generate_incoherent_spectrum_for_angle(0 * EPD, EOR_90, T, Z, 0, eels_photodember, e_w);

[~, ~, psi_incoherent_comb_0] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, EOR_0, T, Z, 0, eels_photodember, e_w);
[~, ~, psi_incoherent_comb_45] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, EOR_0, T, Z, 45, eels_photodember, e_w);
[~, ~, psi_incoherent_comb_90] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, EOR_90, T, Z, 0, eels_photodember, e_w);

%%
setdir = 'article_check/results/';
close all
close all
figure;
image_name = 'tiles_test_meep';
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
close all;
FontName = 'ariel';
FontSize = 15;
clim  = max(abs(EOR_90(:)));
% section_or = 401;
% section_pd = 412;
% create_figure_electricfield(T, Z, EPD, clim/100, setdir, 'field_photodember.png', ...
%     FontSize,section_pd);
create_figure_electricfield(T + 0.2, Z, EOR_90, clim, setdir, 'field_rectification_meep.png', ...
    FontSize);

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

function create_figure_electricfield(T, Z, E, clim, setdir, filename, FontSize, section)
    figure;
    imagesc(T(:,1), Z(1,:), E, [-clim, clim]);
%     h = gcf();
%     h.Renderer = "painters";
    set(gca,'FontSize',FontSize);
    xlim([-.3,1.5]);
    ylim([-100,100]);
    xticks(-.3:.3:1.5)
    colormap(utils.redblue);
    pbaspect([1 1 1])
    set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
    colorbar;
%     l0 = create_line(0); l1 = create_line(0.3); l2 = create_line(0.7); 
%     
    if nargin>7
        hold on
        xx = T(:,1); yyor = Z(:,section);
        plot(xx,yyor,LineStyle="--",Color='#808080',LineWidth=1)
    end

    set(groot,'defaultAxesXTickLabelRotationMode','manual');
    exportgraphics(gcf, [setdir, filename],'resolution', 300);
%     print(h,'-vector', '-dsvg', [setdir, filename]) 
%     saveas(h, [setdir, filename],'resolution', 300) 
    
end
clear all
close all

e_w = linspace(-5,5,181);


[tc, xc, EPD_z] = photodember_meep();
spot_size_fwhm = 80;
[~, ~, field_laser_profile_x, field_laser_profile_y, field_laser_profile_z] = rectification_meep_triple(spot_size_fwhm);
EOR_z = field_laser_profile_z.' * 1e3;
EOR_x = field_laser_profile_x.' * 1e3;

EPD_z = EPD_z * 1.5 - EOR_x / 2;
[TPD, ZPD] = ndgrid(tc, xc);
%%



%%

%%
T = TPD;
Z = ZPD;

for shift = [0]
    EOR = 0.8 * circshift(EOR_z,shift,2);
    EPD = 0.8 * circshift(EPD_z,-shift,2);
    [t_w , psi_incoherent_comb_0] = incoherent_spectrum_from_meep_fields(-EOR + EPD, T, Z, e_w);
    [~ , psi_incoherent_comb_45] = incoherent_spectrum_from_meep_fields(EPD, T, Z, e_w);
    [~ , psi_incoherent_comb_90] = incoherent_spectrum_from_meep_fields(EOR + EPD, T, Z, e_w);
    
    [~ , psi_incoherent_or_0] = incoherent_spectrum_from_meep_fields(-EOR, T, Z, e_w);
    [~ , psi_incoherent_or_45] = incoherent_spectrum_from_meep_fields(-EOR * 0, T, Z, e_w);
    [~ , psi_incoherent_or_90] = incoherent_spectrum_from_meep_fields(EOR, T, Z, e_w);
    
    [~ , psi_incoherent_pd] = incoherent_spectrum_from_meep_fields(EPD, T, Z, e_w);
    
    %%
    
    % eels_photodember = setup_parameters_eels_photodember(10, 40e-6);
    % w = eels_photodember.electron.incoherent_gaussian_blurring_window(e_w,t_w);
    % psi_incoherent =  eels_photodember.incoherent_convolution(psi_assemb, w, t_w, e_w);
    
    
    %% Get measurement data
    
    [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
    [~, ~, ~, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));
    [~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));
    
    
    %%
    
    t_w = t_w - 0.4;
    setdir = 'article_check/results/shift_trials/';
    close all
    close all
    figure;
    image_name = strcat('tiles_check_meep',num2str(shift));
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
end
%%
    
function [t_w , psi_incoherent] = incoherent_spectrum_from_meep_fields(EZ, T, Z, e_w)

    [t_w, eels_spectra] = utils_spectrum.calculate_spectrum_from_fields(EZ, T, Z*1e-6);
    psi_assemb = utils_spectrum.spectrum_to_coherent_eels_mod(e_w, t_w, eels_spectra * 1e2 * 3);
    eels_photodember = setup_parameters_eels_photodember(10, 40e-6);
    w = eels_photodember.electron.incoherent_gaussian_blurring_window(e_w,t_w);
    psi_incoherent =  eels_photodember.incoherent_convolution(psi_assemb, w, t_w, e_w);

end



function [tc, xc, field_pd] = photodember_meep()
    load('article_check/saved_matrices_meep/photodember/spot_size_30_shift04/field_ez_pd_intensity_10t0_0.5.mat')
    

  
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


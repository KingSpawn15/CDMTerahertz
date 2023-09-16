clear all
close all

or_spot_sigma = 50e-6;
pd_spot_fwhm = 40e-6;
pd_z_max = 90e-6;
pulse_energy_experiment_nj = 10;


%% Get measurement    
[mean_exp_0, ~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
[mean_exp_90, ~, ~, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));
[mean_exp_45, deltat, energy, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));

%% Photodember field

[TPD, ZPD, EPD,psi_sub_pd, psi_incoherent_pd, ...
    eels, w, e_w, t_w, tt, zz] = electric_field_photodember(pulse_energy_experiment_nj, ...
    pd_spot_fwhm, pd_z_max);
EPD = EPD*(1.26/6.34);

%% Rectification field
rectification_param_z0 = -1e-6;
rectification_param_thickness = .5e-3;
rectification_param_wavelength = 800e-9;
rectification_param_pulse_time_tau = 30e-15;
rectification_param_pulse_laser_spot_sigma = or_spot_sigma;

[~, ~,TOR, ZOR, EOR] = eels_theoretical_2(rectification_param_pulse_time_tau, ...
    rectification_param_wavelength, ...
    rectification_param_thickness, ...
    rectification_param_z0, ...
    rectification_param_pulse_laser_spot_sigma);

EOR = EOR.';
ZOR = ZOR .* 1e6;

%%
% Define the conditions
condition1 = ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:));
condition2 = TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:));

% Apply conditions to the arrays
EOR1 = EOR(condition1, condition2);
TOR1 = TOR(condition2, condition1);
ZOR1 = ZOR(condition2, condition1);
EORintrap = interp2(TOR1.', ZOR1.', EOR1, TPD.', ZPD.', 'linear', 0);
EPDintrap = interp2(TPD.', ZPD.', EPD, TOR.', ZOR.', 'linear', 0);

%%
% [t0_vec_or, eels_cc_or] = utils_spectrum.calculate_spectrum_from_fields(EOR.', TOR, ZOR* 1e-6);
[t0_vec_or, spectrum_or] = utils_spectrum.calculate_spectrum_from_fields((EOR + EPDintrap).', TOR, ZOR* 1e-6);
%%
% psi_sub_or = utils_spectrum.spectrum_to_coherent_eels(t_w, e_w, spectrum_or, t0_vec_or);
% psi_incoherent_or = eels.incoherent_convolution(psi_sub_or, w, t_w, e_w);
% 
% FontName = 'ariel';
% FontSize = 10;
% figure;
% plot_tile(e_w,t_w, psi_incoherent_or);
% set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
% set(gcf,'position', [200 , 200 , 200 + 150, 200 + 150]);

%%
close all
setdir = 'figure 1 fields/results/';
FontName = 'ariel';
FontSize = 10;

clim = max(abs(EOR1(:)));
create_figure_electricfield(TPD, ZPD, EPD, clim/100, setdir, 'field_photodember.png', FontSize);
create_figure_electricfield(TOR1, ZOR1, EOR1, clim, setdir, 'field_rectification.png', FontSize);
create_figure_electricfield(TPD, ZPD, EORintrap + EPD, clim, setdir, 'field_combined.png', FontSize);

%%

function create_figure_electricfield(T, Z, E, clim, setdir, filename, FontSize)
    figure;
    imagesc(T(:,1), Z(1,:), E, [-clim, clim]);
    set(gca,'FontSize',FontSize);
    xlim([-.3,1.5]);
    ylim([-100,100]);
    xticks(-.3:.3:1.5)
    colormap(utils.redblue);
    pbaspect([2 1 1])
    set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
    colorbar;
    exportgraphics(gcf, [setdir, filename],'resolution', 300);
end



function [TT, ZZ, electric_field_zt,psi_sub_pd ,...
    psi_incoherent_pd,eels, w, e_w, t_w, tt, zz] = electric_field_photodember(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max)

    % set up the eels parameters
    [eels, w, e_w, t_w, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max);
    
    % calculate the photodember interaction voltage
    interact_v_pd_store = utils.eels_pdor(eels,'photodember');
    

    interact_v_pd = circshift(interact_v_pd_store, [0,0]);
    

    alpha_pd_0 =  0.05;
    

    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 ;
    

    [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);

    [TT, ZZ] = ndgrid(tt * 1e12,zz * 1e6);
    

    electric_field_zt = utils_spectrum.derivative_f_dz(interact_v_pd_store, zz, tt);
    electric_field_zt(isnan(electric_field_zt)) = 0;

end

function [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max)

% (TODO: Modify default_parameters)
     [laser_parameters,discretization_params, utem_parameters,...
        numerical_parameters] = default_parameters_2(pd_spot_fwhm);

    laser_parameters.pulse_energy_experiment = 0.1 * pulse_energy_experiment_nj *1e-9;
    discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
    discretization_params.delay_max = 2 * 1.5e-12;
    discretization_params.z_max = pd_z_max;
    
    utem_parameters.electron_total_energy = 0.94;
    laser_parameters.laser_pulse_time_fwhm = 650e-15;
    laser_parameters.theta_pol = 90*pi/180;
    
    laser = Laser(laser_parameters);
    
    discretization = Discretization(discretization_params);
    elec = UTEMElectron(utem_parameters);
    
    [w, e_w, t_w_store] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
        discretization.energy, discretization.deltat);
    
    
    eels_parameters.electron = elec;
    eels_parameters.discretization = discretization;
    eels_parameters.numerical_parameters = numerical_parameters;
    eels_parameters.laser = laser;
    eels_parameters.material = IndiumArsenide();
    eels_parameters.materialmaterial.phase = 0;
    eels = EELS(eels_parameters);

    tt = discretization.t;
    zz = discretization.z;
end

function plot_tile(x, y, z)
    FontSize = 10;
    FontName = 'ariel';
    nexttile
    imagesc(x, y, z);
    ylim([-1,1.5]);
    xlim([-5,5]);
    colormap jet
    axis square
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

% function [mean_eels, deltat, energy, eels_measure] = get_measurement(angle)
% 
% 
%     [eels_measure, energy,...
%         time] = measurement_plots.data_measurement(angle);
%     
%     deltat = time(time > -1 & time < 1.5);
%     eels_measure = eels_measure(time > -1 & time < 1.5,:);
%     
%     x_c = zeros(1, size(deltat,2));
%     
%     eels_measure = eels_measure ./ sqrt(sum(eels_measure.^2, 2));
%     
%     
%     for i = 1:size(eels_measure,1)
%         [~, ind_max] = max(eels_measure(i,:));
%         x_c(i) = ind_max;
%     end
%     
%     mean_eels = energy(floor(x_c));
% end

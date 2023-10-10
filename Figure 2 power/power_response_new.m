clear all
close all

%% Get measurement data
[com,deltat,energy,eels_measure,errs] = getmeasurement_power();
%%
delete(gcp('nocreate'))
% create a new parallel pool with 6 workers
parpool(6);

or_spot_sigma = 50e-6;
pd_spot_fwhm = 40e-6;
pd_z_max = 90e-6;


%Photodember field
pulse_energy_list = [0.1000    0.1700    0.3000    0.5200    1.0000    1.7300    3.0000    5.1900   10.0000];
kfactor = 1.2;
for ii = 1:length(pulse_energy_list)
    disp(ii)
    [TPD, ZPD, EPD{ii},psi_sub_pd{ii}, psi_incoherent_pd{ii}, ...
        eels, w, e_w, t_w, tt, zz] = electric_field_photodember(pulse_energy_list(ii), ...
        pd_spot_fwhm, pd_z_max);
    EPD{ii} = EPD{ii}*(1.26/6.34);
end



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

% Define the conditions
condition1 = ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:));
condition2 = TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:));

% Apply conditions to the arrays
EOR1 = EOR(condition1, condition2);
TOR1 = TOR(condition2, condition1);
ZOR1 = ZOR(condition2, condition1);
% EORintrap = interp2(TOR1.', ZOR1.', EOR1, TPD.', ZPD.', 'linear', 0);
for ii = 1:length(pulse_energy_list)
    EPDintrap{ii} = interp2(TPD.', ZPD.', kfactor * EPD{ii}, TOR.', ZOR.', 'linear', 0);
end

ZOR = ZOR .* 1e-6;
%%
% ii = 1;
for ii = 1:length(pulse_energy_list)
    disp(pulse_energy_list(ii))
    [t0_vec,eels_comb{ii}] = utils_spectrum.calculate_spectrum_from_fields((0.1 * pulse_energy_list(ii)) * EOR + ...
        EPDintrap{ii}, TOR, ZOR);
    psi_incoherent_comb{ii} = eels.incoherent_convolution(utils_spectrum.spectrum_to_coherent_eels(t_w, e_w, eels_comb{ii}, t0_vec),...
        w, t_w, e_w);
end

%%
close all
figure;
FontName = 'ariel';
FontSize = 10;
ttt = tiledlayout(3,3,"TileSpacing","compact");
ttt.Padding = "loose"
nexttile
imagesc(energy, deltat - 0.3, eels_measure{1});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile;
imagesc(energy, deltat- 0.3, eels_measure{8});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile
imagesc(energy, deltat- 0.3, eels_measure{9});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{1});
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
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3])
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{9});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
axis square
nexttile([1,3]);
% errors = 0.1 * ones(size(deltat));
colors = turbo(9);
for i = 1:1:9
    plot(eels_comb{i} + 4 * (i-1), t0_vec, 'color', colors(i,:),'LineWidth',1);
    hold on;
    errorbar(com(i,2:2:end) + 4 * (i-1), deltat(2:2:end) - 0.3, ...
       errs(i,2:2:end), 'horizontal',  'color', colors(i,:),'LineStyle','none','LineWidth',.5);
set(gca, 'YDir', 'reverse');
    
end
hold off
set_axis_properties(gca,FontSize,FontName,.2,-1:0.2:1.5,0:4:32,'','',FontSize,[0 0 0])
% axis equal;
% daspect([1 1 1]);
ylim([-0.5, 0.5]);
% xlim([-1,10]);

set(gcf,'Position',[200,200,200 + 450,200 + 600]); %set paper size (does not affect display)

exportgraphics(gcf, 'figure 2 power/fig2_new.png', 'Resolution',300);


%%

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

function psi_incoherent = generate_incoherent_spectrum_for_angle(EPD, EOR, TOR, ZOR, eels, w, t_w, e_w, theta)
    
    [t0_vec,eels_comb] = utils_spectrum.calculate_spectrum_from_fields(-cos(2 * theta * pi / 180) * EOR + ...
        EPD, TOR, ZOR);
    psi_incoherent =  eels.incoherent_convolution(utils_spectrum.spectrum_to_coherent_eels(t_w, e_w, eels_comb, t0_vec),...
        w, t_w, e_w);

end


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

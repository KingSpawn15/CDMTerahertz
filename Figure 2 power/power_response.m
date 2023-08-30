clear all;
close all;
delete(gcp('nocreate'))
delete(gcp('nocreate'))
% parpool(6);




%%
pulse_energy_list = [0.1000    0.1700    0.3000    0.5200    1.0000    1.7300    3.0000    5.1900   10.0000];

i = 1;
for pulse_energy = pulse_energy_list
    [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy);
    interact_v_pd_store_cell{i} = utils.eels_pdor(eels,'photodember');
    i = i+1;
end

interact_v_or_store = utils.eels_pdor(eels,'rectification');


%%

i = 1;
for pulse_energy = pulse_energy_list
    alpha_pd_0 =  .07;
    alpha_or_0 = 20 * 1.3 * 1.2 * pulse_energy / 10;
    
    interact_v_pd = circshift(interact_v_pd_store_cell{i}, [18,0]);
    interact_v_or = circshift(interact_v_or_store, [-15,0]);
    
    t_w = t_w_store-0.2;
    
    
    % alpha_pd = alpha_pd_0; 
    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 + interact_v_or * alpha_or_0;
    [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);
    
    % % alpha_pd = alpha_pd_0; 
    % loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 + interact_v_or * alpha_or_0;
    % [psi_sub_comb , psi_incoherent_comb] = eels.energy_loss_spectrum(loss_spectrum_parameters);
    

    [~, eels_ind_pd] = max(psi_sub_pd,[],2);
    e_exc_pd{i} = e_w(eels_ind_pd);
    psi_incoherent_pd_cell{i} = psi_incoherent_pd;
    i = i+1;
end



[com,deltat,energy,eels_measure] = getmeasurement2();



%%
close all
figure;
FontName = 'ariel';
FontSize = 10;
ttt = tiledlayout(3,3,"TileSpacing","compact");
ttt.Padding = "loose"
nexttile
imagesc(energy, deltat-0.3, eels_measure{1});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile;
imagesc(energy, deltat-0.3, eels_measure{8});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile
imagesc(energy, deltat-0.3, eels_measure{9});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile
imagesc(e_w,t_w, psi_incoherent_pd_cell{1});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3])
nexttile
imagesc(e_w,t_w, psi_incoherent_pd_cell{8});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3])
nexttile
imagesc(e_w,t_w, psi_incoherent_pd_cell{9});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
axis square
nexttile([1,3]);
errors = 0.1 * ones(size(deltat));
colors = turbo(9);
for i = 1:1:9
    plot(0.25 * e_exc_pd{i} + 1 * (i-1), t_w, 'color', colors(i,:),'LineWidth',1);
    hold on;
    errorbar(0.25 * com(i,2:1:end) + 1 * (i-1), deltat(2:1:end) - 0.3, ...
        .5 * errors(2:1:end), 'horizontal',  'color', colors(i,:),'LineStyle','none','LineWidth',1);
set(gca, 'YDir', 'reverse');
    
end
hold off
set_axis_properties(gca,FontSize,FontName,.2,-1:0.2:1.5,0:1:8,'','',FontSize,[0 0 0])
% axis equal;
% daspect([1 1 1]);
ylim([-0.5, 0.5]);
xlim([-1,10]);

set(gcf,'Position',[200,200,200 + 450,200 + 600]); %set paper size (does not affect display)

exportgraphics(gcf, 'figure 2 power/fig2.png', 'Resolution',300);



%% Function to set axis properties

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end

function [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj)
%     spot_fwhm = 40e-6;
%     [laser_parameters,discretization_params, utem_parameters,...
%         numerical_parameters] = default_parameters_2(spot_fwhm);

    
    if nargin < 1
        pulse_energy_experiment_nj = 10;
    end

     [laser_parameters,discretization_params, utem_parameters,...
        numerical_parameters] = default_parameters_2();

    laser_parameters.pulse_energy_experiment = 0.1 * pulse_energy_experiment_nj *1e-9;
    discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
    discretization_params.delay_max = 2 * 1.5e-12;
    discretization_params.z_max = 30e-6;
    
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



function [com,deltat, energy, eels_measure] = getmeasurement2()
    load('saved_matrices\PulseEnergy.mat')
    deltat = Time;
    energy = EnergyCrop;
    % assign variables deltat and energy to values of Time and EnergyCrop, respectively
    % these variables will be used later in the code
    
   
    
    x_c = zeros(11, size(Time,2));
    
    for ii = 1:11
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        eels_measure{ii} = eels;
        % Define the mass of each element in the matrix (assuming unit mass)
        % mass_matrix = ones(size(eels));
        mass_matrix = eels;
        % Initialize a vector to store the x coordinates of the center of mass of each row
       
    
        % Calculate the x coordinate of the center of mass of each row
%         for i = 1:size(eels,1)
%             row_mass = sum(mass_matrix(i,:));
%             for j = 1:size(eels,2)
%                 x_c(ii,i) = x_c(ii,i) + j*mass_matrix(i,j);
%             end
%             x_c(ii,i) = x_c(ii,i)/row_mass;
%         end

        for i = 1:size(eels,1)
%             row_mass = sum(mass_matrix(i,:));
            for j = 1:size(eels,2)
                [~, ind_max] = max(mass_matrix(i,:));
                x_c(ii,i) = ind_max;
            end
            x_c(ii,i) = x_c(ii,i);
        end
    end

    com = energy(floor(x_c));
end

clearvars;
dir = 'results/optimization_results/';
[~ , ~ , ~] = mkdir(dir);
load('saved_matrices/v_struct_5.mat');
load('results/optimization/optimization_results.mat');

[status, msg, msgID] = mkdir('results/optimization_plots');

% best parameter 1
% x = [0.0109 0.4259 0]
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters();

discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;
discretization_params.delay_max = 1.5e-12;

laser = Laser(laser_parameters);
discretization = Discretization(discretization_params);
% elec = UTEMElectron(utem_parameters);
%
% [w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
%     discretization.energy, discretization.deltat);
% [fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);

% eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.material = IndiumArsenide();
eels_parameters.numerical_parameters = numerical_parameters;

nplots = 20;
list_weights_pd = [optimization_combined(1:nplots).weight_pd];
list_weights_or = [optimization_combined(1:nplots).weight_or];
list_delay = [optimization_combined(1:nplots).delay];
list_electron_energy = [optimization_combined(1:nplots).electron_tot_energy];
list_electron_theta = [optimization_combined(1:nplots).electron_theta_rad];

for ii = 1:nplots
    
    utem_parameters.electron_total_energy = list_electron_energy(ii);
    utem_parameters.electron_total_time_fs = 350;
    utem_parameters.electron_time_coherent_fwhm_fs = 20;
    utem_parameters.electron_theta = list_electron_theta(ii);
    
    delay = list_delay(ii);
    interaction_gain_factor_photodember = list_weights_pd(ii);
    interaction_gain_factor_rectification = list_weights_or(ii);
    
    plot_ind = 1;
    close all;
    figure = tiledlayout(2,9,'Padding', 'none', 'TileSpacing', 'compact');
    kk = 1;
    
    for theta_pol_degree = 10 : 10 : 180
        
        
        elec = UTEMElectron(utem_parameters);
        laser.theta_pol =  theta_pol_degree.*(pi/180);
        eels_parameters.laser = laser;
        eels_parameters.electron = elec;
        
        eels = EELS(eels_parameters);
        
        [w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
            discretization.energy, discretization.deltat);
        
        loss_spectrum_parameters.method = 'combination';
        loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
            v_struct.(strcat('angle_',num2str(theta_pol_degree))) + ...
            interaction_gain_factor_photodember * circshift(v_struct.(strcat('photodember')),[delay 0]);
        tic;
        [psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);
        toc;
        nexttile;
        imagesc(e_w,t_w, psi_incoherent);
        
        ylim([-1 , 1.5]);

        if plot_ind == 9 || plot_ind == 18
            colorbar;
        end
        colormap jet
        axis square
        ax = gca;
        ax.FontSize = 10;
        ax.LineWidth = 1;
        if (plot_ind == 1 || plot_ind == 10)
            ax.YTick = -1:0.5:1.5;
            ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',12);
        else
            ax.YTick = [];
        end
        
        if (plot_ind >= 10)
            ax.XTick = -4:2:4;
            xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',12);
        else
            ax.XTick = [];
        end
        
        box on
        plot_ind = plot_ind + 1;
        
    end
    set(gcf,'Position',[100,100,1600,400]);
    str = [dir,'no_',num2str(ii),'_',...
        'pd_gain=',num2str(interaction_gain_factor_photodember,'%.2f'),...
        'or_gain=',num2str(interaction_gain_factor_rectification,'%.2f'),...
        'delay=',num2str(delay),...
        'elec_eng=',num2str(utem_parameters.electron_total_energy,'%.2f'),...
        'elec_theta=',num2str(utem_parameters.electron_theta*180/pi,'%.2f')
        ];
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
end
close all
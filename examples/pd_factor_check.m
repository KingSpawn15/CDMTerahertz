clearvars;

for energy_gain_factor = 1
    %     energy_gain_factor = 1;
    dir = 'results/or_factor_check/';
    [~ , ~ , ~] = mkdir(dir);
    load('saved_matrices/v_struct_5deg.mat');
    factor_rect = 0.3;
    factor_pd = 0;
    delay_pd = -15;
    delay_or = 6;
    
    [laser_parameters,discretization_params, utem_parameters,...
        numerical_parameters] = default_parameters_2();
    
    discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;
    discretization_params.delay_max = 1.5e-12;
    
    utem_parameters.electron_total_energy = 0.94;
    utem_parameters.electron_total_time_fs = 350;
    utem_parameters.electron_time_coherent_fwhm_fs = 20;
    utem_parameters.electron_theta = -6.35*pi/180;
    utem_parameters.electron_velocity_c = 0.7;
    numerical_parameters.subsampling_factor = 20;
    
    pulse_width_fs = 50;
    % energy_gain_factor = 50;
    laser_parameters.pulse_energy_gain_factor = 0.014 * energy_gain_factor;
    laser_parameters.laser_pulse_time_fwhm = pulse_width_fs*1e-15;
    laser = Laser(laser_parameters);
    
    discretization = Discretization(discretization_params);
    elec = UTEMElectron(utem_parameters);
    
    [w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
        discretization.energy, discretization.deltat);
    [fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);
    
    eels_parameters.electron = elec;
    eels_parameters.discretization = discretization;
    material_n = IndiumArsenide();
    % material_n.phase = 0;
    % material_n.gamma = 3.51e12;%[s^-1]
    % material_n.gamma_factor = 1;
    eels_parameters.material = material_n;
    
    eels_parameters.numerical_parameters = numerical_parameters;
    
    
    
    plot_ind = 1;
    close all;
    figure = tiledlayout(1,1,'Padding', 'none', 'TileSpacing', 'compact');
    kk = 1;
    theta_pol_degree = 90;
    
    
    laser.theta_pol =  theta_pol_degree.*(pi/180);
    eels_parameters.laser = laser;
    
    eels = EELS(eels_parameters);
    
    loss_spectrum_parameters.interaction_gain_factor_rectification = ...
        1;
    loss_spectrum_parameters.interaction_gain_factor_photodember =...
        0 ;
    loss_spectrum_parameters.method = 'rectification';
    v_or = eels.interaction_v(loss_spectrum_parameters);
    
    loss_spectrum_parameters.interaction_gain_factor_rectification = ...
        0;
    loss_spectrum_parameters.interaction_gain_factor_photodember =...
        1 ;
    loss_spectrum_parameters.method = 'photodember';
    v_pd = eels.interaction_v(loss_spectrum_parameters);
    
    
    for interaction_gain_factor_rectification = [0.2, 0.5, 1, 1.5, 2]
        for interaction_gain_factor_photodember = factor_pd 
            
            loss_spectrum_parameters.interaction_gain_factor_rectification = ...
                interaction_gain_factor_rectification;
            loss_spectrum_parameters.interaction_gain_factor_photodember =...
                factor_pd * interaction_gain_factor_photodember ;
            loss_spectrum_parameters.method = 'rectification';
            %             loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
            %                 circshift(v_struct.(strcat('angle_',num2str(theta_pol_degree))),[0  0]) + ...
            %                 interaction_gain_factor_photodember * circshift(v_struct.(strcat('photodember')),[delay 0]);
            loss_spectrum_parameters.interact_v = factor_rect * interaction_gain_factor_rectification * ...
                circshift(v_or,[delay_or  0]) + ...
                interaction_gain_factor_photodember * circshift(v_pd,[delay_pd 0]);
            
            [psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);
%             nexttile;
            imagesc(e_w,t_w, psi_incoherent);
            
            ylim([-1 , 1.5]);
            
            %                 xlim([-3,3]);
            %                 ylim([-0.5,1]);
            
            %                 if plot_ind == 9 || plot_ind == 18
            if 1==1
                colorbar;
            end
            colormap jet
            axis square
            ax = gca;
            ax.FontSize = 18;
            ax.LineWidth = 1;
            %             if (plot_ind == 1 || plot_ind == 10)
            if 1==1
                ax.YTick = -1:0.5:1.5;
                ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
            else
                ax.YTick = [];
            end
            
            %             if (plot_ind >= 10)
            if 1==1
                ax.XTick = -4:2:4;
                xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
            else
                ax.XTick = [];
            end
            
            box on
            plot_ind = plot_ind + 1;
            
            str = [dir,'combination_',...
                'energy_gain=',num2str(energy_gain_factor),...
                'angle=',num2str(theta_pol_degree),...
                'pd_gain=',num2str(interaction_gain_factor_photodember),...
                'or_gain=',num2str(interaction_gain_factor_rectification)
                ];
            exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
            
        end
        
        
    end
end






clearvars;
[~ , ~ , ~] = mkdir('results/combination/exhaustive/');
load('saved_matrices/v_struct_5.mat');
[status, msg, msgID] = mkdir('results');

% best parameter 1
% x = [0.0109 0.4259 0]
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters();

discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;
discretization_params.delay_max = 1.5e-12;
discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;
discretization_params.delay_max = 1.5e-12;

utem_parameters.electron_total_energy = 0.96;
utem_parameters.electron_total_time_fs = 120;
utem_parameters.electron_time_coherent_fwhm_fs = 20;
utem_parameters.electron_theta = -3*pi/180;
utem_parameters.electron_velocity_c = 0.7;


% utem_parameters.electron_total_energy = 0.8;
% utem_parameters.electron_total_time_fs = 150;
% utem_parameters.electron_time_coherent_fwhm_fs = 20;
% utem_parameters.electron_theta = -5*pi/180;
% utem_parameters.electron_velocity_c = 0.7;


laser = Laser(laser_parameters);
discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);
[fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);

eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.material = IndiumArsenide();
eels_parameters.numerical_parameters = numerical_parameters;

% load('+optimization/last_optimal_xarr.mat','xarr');
% factor_pd = xarr(3).x(1);
% factor_or = xarr(3).x(2);
% delay = xarr(3).x(3);
% 
% factor_pd = -.17;
% factor_or = .07; 
% delay = 10;
for delay = 10:5:60
for interaction_gain_factor_rectification = -0.1:0.01:0.15
    for interaction_gain_factor_photodember =-0.1:0.01:0.15
        phase = 0;
        plot_ind = 1;
        close all;
        figure = tiledlayout(2,9,'Padding', 'none', 'TileSpacing', 'compact');
        kk = 1;
        for theta_pol_degree = 10:10:180
            
            laser.theta_pol =  theta_pol_degree.*(pi/180);
            eels_parameters.laser = laser;
            
            eels = EELS(eels_parameters);
            
            
            loss_spectrum_parameters.interaction_gain_factor_rectification = ...
                interaction_gain_factor_rectification;
            loss_spectrum_parameters.interaction_gain_factor_photodember =...
                interaction_gain_factor_photodember * exp(1i*phase);
            loss_spectrum_parameters.method = 'combination';
            loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
                v_struct.(strcat('angle_',num2str(theta_pol_degree))) + ...
                interaction_gain_factor_photodember * circshift(v_struct.(strcat('photodember')),[delay 0]);
            [psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);
            
            nexttile;
            imagesc(e_w,t_w, psi_incoherent);

            ylim([-1 , 1.5])
            
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
        str = ['results/combination/exhaustive/','combination_',...
            'pd_gain=',num2str(interaction_gain_factor_photodember),...
            'or_gain=',num2str(interaction_gain_factor_rectification),...
            'delay=',num2str(delay)];
                exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
        
    end
end
end




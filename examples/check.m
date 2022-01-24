clearvars;
[~ , ~ , ~] = mkdir('results/combination');
% load('saved_matrices/v_struct_3.mat');
[status, msg, msgID] = mkdir('results');

% best parameter 1
% x = [0.0109 0.4259 0]
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters();
discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;


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

lst = [[0];[0]];
for interaction_gain_factor_rectification = [0.1]
    for interaction_gain_factor_photodember = 0
        phase = 0;
        plot_ind = 1;
        close all;
        figure = tiledlayout(2,1,'Padding', 'none', 'TileSpacing', 'compact');
        
        for iter = lst
            
            theta_pol_degree = iter(1);
            delay = iter(2);
            
            laser.theta_pol =  theta_pol_degree.*(pi/180)*0;
            eels_parameters.laser = laser;
            
            eels = EELS(eels_parameters);
            
            
            loss_spectrum_parameters.interaction_gain_factor_rectification = ...
                interaction_gain_factor_rectification;
            loss_spectrum_parameters.interaction_gain_factor_photodember =...
                interaction_gain_factor_photodember * exp(1i*phase);
            loss_spectrum_parameters.method = 'rectification';
            %             loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
            %                 v_struct.(strcat('angle_',num2str(theta_pol_degree))) + ...
            %                 exp(1i*phase) * interaction_gain_factor_photodember * v_struct.(strcat('photodember'));
            
%             loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
%                 v_struct.(strcat('angle_',num2str(theta_pol_degree))) + ...
%                 interaction_gain_factor_photodember * circshift(v_struct.(strcat('angle_',num2str(theta_pol_degree))),[delay 0]);
%             
            %             loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
            %                 v_struct.(strcat('angle_',num2str(theta_pol_degree)));
            
            
            [psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);
            
            nexttile;
            imagesc(e_w,t_w, psi_incoherent);
            
            ylim([-1 , 1.5])
            
            if plot_ind == 6 || plot_ind == 12
                colorbar;
            end
            colormap jet
            axis square
            ax = gca;
            ax.FontSize = 14;
            ax.LineWidth = 1;
            if (plot_ind == 1 || plot_ind == 7)
                ax.YTick = -1:0.5:2.5;
                ylabel('\Deltat [ps]')
            else
                ax.YTick = [];
            end
            
            if (plot_ind >= 7)
                xlabel('Energy [eV]')
            else
                %                 ax.XTick = [];
            end
            
            box on
            plot_ind = plot_ind + 1;
        end
        set(gcf,'Position',[100,100,1200,400]);
        str = ['results/combination/','combination',...
            'pd_gain=',num2str(interaction_gain_factor_photodember),...
            'or_gain=',num2str(interaction_gain_factor_rectification)];
        %         exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
        
    end
end


clear vars;
load('saved_matrices/v_struct.mat');

[status, msg, msgID] = mkdir('results');


[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters();

laser = Laser(laser_parameters);
discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);


interaction_gain_factor = 1e-1;
interaction_gain_factor_photodember = 0;
 
[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);
[fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);

eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.material = IndiumArsenide();
eels_parameters.numerical_parameters = numerical_parameters;

for interaction_gain_factor_photodember = [0]
    for method = [ "rectification"]
        for theta_pol_degree = 0
            
            laser.theta_pol =  theta_pol_degree.*(pi/180);
            eels_parameters.laser = laser;
            
            eels = EELS(eels_parameters);
            
            loss_spectrum_parameters.method = method;
            loss_spectrum_parameters.interaction_gain_factor = ...
                interaction_gain_factor;
            loss_spectrum_parameters.interaction_gain_factor_photodember =...
                interaction_gain_factor_photodember;
            loss_spectrum_parameters.interact_v = interaction_gain_factor * ...
                v_struct.(strcat('angle_',num2str(theta_pol_degree)));
            
            [psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);
            
            figure(3)
            imagesc(e_w,t_w,psi_sub)
            xlabel('Energy [eV]')
            ylabel('\Deltat [ps]')
            colorbar
            colormap jet
            axis square
            drawnow
            
            close all;
            figure(4)
            imagesc(e_w,t_w, psi_incoherent)
            xlabel('Energy [eV]')
            ylabel('\Deltat [ps]')
            ylim([-1 , 1.5])
            colorbar
            colormap jet
            axis square
            ax = gca;
            ax.FontSize = 18;
            ax.LineWidth = 1;
            ax.YTick = -1:0.5:2.5;
            box on
            drawnow
            
            str = ['results/','eels_angle=',num2str(theta_pol_degree),'_method=',char(method),...
                'pd_gain=',num2str(interaction_gain_factor_photodember)];
            %             savefig(gcf,strcat(str,'.fig'));
%             exportgraphics(gcf,strcat(str,'.png'),'Resolution',300);
        end
    end
end
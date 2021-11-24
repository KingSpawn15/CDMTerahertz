% for angle = [90:4:180]
%     close all
%     measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
%     str = ['results/measurement_plots/','eels_angle=',num2str(angle)];
%     savefig(gcf,strcat(str,'.fig'));
%     exportgraphics(gcf,strcat(str,'.png'),'Resolution',300);
% end
%


% electron parameters

electron_total_energy = 0.8;
electron_total_time_fs = 150;
electron_time_coherent_fwhm_fs = 20;
electron_theta = -5*pi/180;
electron_velocity_c = 0.7;

% laser parameters
pulse_energy_experiment = 10 * 1e-9;
pulse_energy_gain_factor = 0.05;
% theta_pol_degree = 90;
% theta_pol = theta_pol_degree.*(pi/180);

% Discretization

ddt = 10e-15;   delay_max = 2.5e-12;
ddz = 1e-6;     zmax = 1e-4;
fs = 2.4e15;    l = 2.4e4;

% subsampling
subsampling_factor = 60;

% sample geometrical parameters
x0 = 0;
y0 = -1e-6;

% effect
% method = 'rectification';

% fitting parameters
interaction_gain_factor = 1e-1;

elec = UTEMElectron(electron_total_energy,...
    electron_total_time_fs, electron_time_coherent_fwhm_fs,...
    electron_theta, electron_velocity_c);

discretization = DiscretizationSteps(ddt, delay_max, ddz, zmax, fs, l);

sample_parameters = Sample(x0, y0);

for interaction_gain_factor_photodember = [0.01, 0.01, 1]
    for method = [ "combination" , "rectification" ]
        for theta_pol_degree = [0:5:90]
            theta_pol = theta_pol_degree.*(pi/180);
            las = Laser(pulse_energy_experiment, pulse_energy_gain_factor, theta_pol);
            eels = EELS(elec, las, discretization, sample_parameters);
            
            
            [w, e_w, t_w] = elec.subsampling(subsampling_factor,...
                discretization.energy, discretization.deltat);
            [fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);
            
            % calculate dynamics
            
            interact_v = eels.interaction_v(method, interaction_gain_factor,...
                interaction_gain_factor_photodember);
            f_t = eels.calculate_ft(interact_v);
            psi_coherent = eels.calculate_psi_coherent(f_t);
            psi_sub = EELS.psi_sub_sampled(subsampling_factor, psi_coherent , e_w);
            
            figure(3)
            imagesc(e_w,t_w,psi_sub)
            xlabel('Energy [eV]')
            ylabel('\Deltat [ps]')
            colorbar
            colormap jet
            axis square
            drawnow
            
            psi_incoherent = EELS.incoherent_convolution(psi_sub, w, t_w, e_w);
            
            close all;
            figure(4)
            imagesc(e_w,t_w, psi_incoherent)
            xlabel('Energy [eV]')
            ylabel('\Deltat [ps]')
            colorbar
            colormap jet
            axis square
            drawnow
            
            str = ['results/','eels_angle=',num2str(theta_pol_degree),'_method=',char(method),...
                'pd_gain=',num2str(interaction_gain_factor_photodember)];
            savefig(gcf,strcat(str,'.fig'));
            exportgraphics(gcf,strcat(str,'.png'),'Resolution',300);
            
        end
    end
end

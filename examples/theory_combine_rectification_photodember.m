clear all;
laser_parameters.pulse_energy_experiment = 10 * 1e-9;
laser_parameters.pulse_energy_gain_factor = 0.05;
laser_parameters.laser_spot_fwhm = 40e-6;
laser_parameters.theta_pol = 45*pi/180;
laser_parameters.laser_pulse_time_fwhm = 50e-15;
laser = Laser(laser_parameters);

discretization_params.x0 = 0;
discretization_params.y0 = -1e-6;
discretization_params.ddt = 10e-15;   discretization_params.delay_max = 2.5e-12;
discretization_params.fs = 2.4e15;    discretization_params.l = 2.4e4;
discretization_params.t0 = -0.5e-12;
discretization_params.xprime_max = round(3 * laser.laser_spot_sigma,5);
discretization_params.d_xprime = 5e-2 * round(3 * laser.laser_spot_sigma,5);
discretization_params.yprime_max = 1e-6;
discretization_params.d_yprime = 5e-2 * 1e-6;
discretization_params.zprime_max = discretization_params.xprime_max;
discretization_params.d_zprime = discretization_params.d_xprime;
discretization_params.ddz = 1e-6;     discretization_params.zmax = 1e-4;
discretization_params.z_max = 30e-6;
discretization = Discretization(discretization_params);

utem_parameters.electron_total_energy = 0.8;
utem_parameters.electron_total_time_fs = 150;
utem_parameters.electron_time_coherent_fwhm_fs = 20;
utem_parameters.electron_theta = -5*pi/180;
utem_parameters.electron_velocity_c = 0.7;
elec = UTEMElectron(utem_parameters);

% % subsampling
subsampling.tc_subsampling = 30;
subsampling.subsampling_factor = 60;

% material
material = IndiumArsenide();

% fitting
interaction_gain_factor = 1e-1;
% interaction_gain_factor_photodember = 0.5;
% method = "rectification";



[w, e_w, t_w] = elec.subsampling(subsampling.subsampling_factor,...
    discretization.energy, discretization.deltat);
[fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);

for interaction_gain_factor_photodember = [0 , 0.2, 0.5, 1]
    for method = [ "combination" , "rectification" ]
        for theta_pol_degree = 0:15:180
            
            laser.theta_pol =  theta_pol_degree.*(pi/180);
            eels = EELS(elec, laser, discretization , material , subsampling);
            
            interact_v = eels.interaction_v(method, interaction_gain_factor,...
                interaction_gain_factor_photodember);
            
            f_t = eels.calculate_ft(interact_v);
            psi_coherent = eels.calculate_psi_coherent(f_t);
            psi_sub = EELS.psi_sub_sampled(subsampling.subsampling_factor, psi_coherent , e_w);
            
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
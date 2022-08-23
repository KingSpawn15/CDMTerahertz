clearvars;

%%
export_dir = 'examples/phase_study/results/';
base_filename = 'eels_';
pol_angle = 150;
export_file_name = strcat(export_dir,base_filename,num2str(pol_angle),'.png');

%%
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();

discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;
discretization_params.delay_max = 1.5e-12;

utem_parameters.electron_total_energy = 0.94;

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);


eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
laser_parameters.laser_pulse_time_fwhm = 500e-15;
% laser_parameters.theta_pol = pol_angle * pi / 180;
laser_parameters.pulse_energy_experiment = 1e-12;
laser = Laser(laser_parameters);
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
% eels_parameters.material.phase = 0;
% eels_parameters.material.gamma_factor = 0.8;
eels = EELS(eels_parameters);
interaction_gain_factor_rectification = 1;
for interaction_gain_factor_photodember = [0]
    for method =  "rectification"
        for theta_pol_degree = [0 45 90 135]
            
            laser.theta_pol =  theta_pol_degree.*(pi/180);
            eels_parameters.laser = laser;
            eels = EELS(eels_parameters);
            
            loss_spectrum_parameters.method = method;
            loss_spectrum_parameters.interaction_gain_factor_rectification = ...
                interaction_gain_factor_rectification;
            loss_spectrum_parameters.interaction_gain_factor_photodember =...
                interaction_gain_factor_photodember;
            
            v_struct.(strcat('angle_',num2str(theta_pol_degree))) = ...
                eels.interaction_v(loss_spectrum_parameters);
        end
    end
end

method = "photodember";
interaction_gain_factor_photodember = 1;
loss_spectrum_parameters.method = method;
loss_spectrum_parameters.interaction_gain_factor_rectification = ...
    interaction_gain_factor_rectification;
loss_spectrum_parameters.interaction_gain_factor_photodember =...
    interaction_gain_factor_photodember;
v_struct.(strcat('photodember')) = ...
                eels.interaction_v(loss_spectrum_parameters);
%%
for interaction_gain_factor_photodember = [0]
    for method =  "rectification"
        for theta_pol_degree = [150]
            
            laser.theta_pol =  theta_pol_degree.*(pi/180);
            eels_parameters.laser = laser;
            eels = EELS(eels_parameters);
            
            loss_spectrum_parameters.method = method;
            loss_spectrum_parameters.interaction_gain_factor_rectification = ...
                1;
            loss_spectrum_parameters.interaction_gain_factor_photodember =...
                0;
            
            v_struct.(strcat('angle_',num2str(theta_pol_degree))) = ...
                eels.interaction_v(loss_spectrum_parameters);
        end
    end
end

save('examples/phase_study/saved_matrices/v_struct_pdmod.mat','v_struct');

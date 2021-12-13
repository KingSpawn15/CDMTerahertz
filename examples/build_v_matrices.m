clear vars;
[status, msg, msgID] = mkdir('results');

[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters();

laser = Laser(laser_parameters);
discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);


% fitting
interaction_gain_factor = 1;

% [w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
%     discretization.energy, discretization.deltat);
% [fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);

eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.material = IndiumArsenide();
eels_parameters.numerical_parameters = numerical_parameters;

for interaction_gain_factor_photodember = [0]
    for method =  "rectification"
        for theta_pol_degree = 0:15:180
            
            laser.theta_pol =  theta_pol_degree.*(pi/180);
            eels_parameters.laser = laser;
            eels = EELS(eels_parameters);
            
            loss_spectrum_parameters.method = method;
            loss_spectrum_parameters.interaction_gain_factor = ...
                interaction_gain_factor;
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
loss_spectrum_parameters.interaction_gain_factor = ...
    interaction_gain_factor;
loss_spectrum_parameters.interaction_gain_factor_photodember =...
    interaction_gain_factor_photodember;
v_struct.(strcat('photodember')) = ...
                eels.interaction_v(loss_spectrum_parameters);
            
save('saved_matrices/v_struct.mat','v_struct');

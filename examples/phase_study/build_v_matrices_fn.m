function build_v_matrices_fn(user_filename,params)

%%
export_dir = 'examples/phase_study/saved_matrices/';
[~,~,~] = mkdir(export_dir);
base_filename = 'v_struct_';
export_filename = strcat(export_dir, base_filename, user_filename,'.mat');
%%


laser_parameters = params.laser_parameters;
discretization_params = params.discretization_params;
utem_parameters = params.utem_parameters;
numerical_parameters = params.numerical_parameters;

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);
laser = Laser(laser_parameters);

eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
eels = EELS(eels_parameters);

interaction_gain_factor_rectification = 1;
for interaction_gain_factor_photodember = [0]
    for method =  "rectification"
        for theta_pol_degree = [[0: 10 : 180],[45,135]]

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


save(export_filename,'v_struct');

end
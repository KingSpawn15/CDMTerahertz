clearvars;

[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();

numerical_parameters.tc_subsampling = 30;

params.laser_parameters = laser_parameters;
params.discretization_params = discretization_params;
params.utem_parameters =  utem_parameters;
params.numerical_parameters = numerical_parameters;

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);
material = IndiumArsenide();
laser = Laser(laser_parameters);

tic;
interaction_v =  ChargeDynamics.interaction_potential_photodember(discretization, material,...
    laser , numerical_parameters);
toc;
% interact_v_fft = fftshift(fft(interaction_v,length(discretization.t),2),2);
%             interact_v_fft = interact_v_fft./max(discretization.omega);
% 
% eels_parameters.numerical_parameters = numerical_parameters;
% eels_parameters.electron = elec;
% eels_parameters.discretization = discretization;
% eels_parameters.laser = laser;
% eels_parameters.material = IndiumArsenide();
% eels = EELS(eels_parameters);
% % f_t = eels.calculate_ft(interaction_v);
% 
% method = "photodember";
% interaction_gain_factor_photodember = 1e-1;
% loss_spectrum_parameters.method = method;
% 
% loss_spectrum_parameters.interaction_gain_factor_rectification = 0;
% loss_spectrum_parameters.interaction_gain_factor_photodember =...
%     interaction_gain_factor_photodember;
% 
% [psi_sub_com , psi_incoherent_com] = eels.energy_loss_spectrum(loss_spectrum_parameters);
% 
% [w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
%         discretization.energy, discretization.deltat);
% imagesc(e_w,t_w, psi_sub_com);
% % ylim([-1 , 1.5]);


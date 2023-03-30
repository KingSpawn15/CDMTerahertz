clearvars;

[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();

laser_parameters.pulse_energy_experiment = 1e-9;
discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
discretization_params.delay_max = 2 * 1.5e-12;

utem_parameters.electron_total_energy = 0.94;
laser_parameters.laser_pulse_time_fwhm = 650e-15;
laser_parameters.theta_pol = 90*pi/180;

params.laser_parameters = laser_parameters;
params.discretization_params = discretization_params;
params.utem_parameters =  utem_parameters;
params.numerical_parameters = numerical_parameters;

% numerical_parameters.tc_subsampling = 30;

params.laser_parameters = laser_parameters;
params.discretization_params = discretization_params;
params.utem_parameters =  utem_parameters;
params.numerical_parameters = numerical_parameters;

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);
material = IndiumArsenide();
laser = Laser(laser_parameters);

tic;
% interaction_v =  ChargeDynamics.interaction_potential_photodember(discretization, material,...
%     laser , numerical_parameters);

% nteraction_potential_rectification(discretization, material,...
%                 laser , electron, numerical_parameters)

interaction_v =  ChargeDynamics_test.interaction_potential_rectification(discretization, material,...
    laser , elec, numerical_parameters);

toc;
% interact_v_fft = fftshift(fft(interaction_v,length(discretization.t),2),2);
%             interact_v_fft = interact_v_fft./max(discretization.omega);

eels_parameters.numerical_parameters = numerical_parameters;
eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();

eels = EELS(eels_parameters);
% % f_t = eels.calculate_ft(interaction_v);
% 
method = "photodember";
interaction_gain_factor_photodember = 31.2;
loss_spectrum_parameters.method = method;

loss_spectrum_parameters.interaction_gain_factor_rectification = 1;
loss_spectrum_parameters.interaction_gain_factor_photodember =...
    interaction_gain_factor_photodember;
loss_spectrum_parameters.interact_v = interaction_v * interaction_gain_factor_photodember;

[psi_sub_com , psi_incoherent_com] = eels.energy_loss_spectrum(loss_spectrum_parameters);

% 
[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
        discretization.energy, discretization.deltat);
imagesc(e_w,t_w, psi_incoherent_com);
ylim([-1 , 1.5]);


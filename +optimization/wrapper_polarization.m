function [psi_incoherent, e_w, t_w] = wrapper_polarization(params)
%WRAPPER Summary of this function goes here
%   Detailed explanation goes here

interaction_gain_factor_photodember = params.interaction_gain_factor_photodember;
interaction_gain_factor_rectification = params.interaction_gain_factor_rectification;
delay = floor(params.delay);

theta_pol_degree = params.theta_pol_degree;

[~ , ~ , ~] = mkdir('results/combination');
load('saved_matrices/v_struct.mat');
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = optimization.default_parameters();

laser = Laser(laser_parameters);
discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

[ ~ , e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);


eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.material = IndiumArsenide();
eels_parameters.numerical_parameters = numerical_parameters;

laser.theta_pol =  theta_pol_degree.*(pi/180);
eels_parameters.laser = laser;

eels = EELS(eels_parameters);


loss_spectrum_parameters.interaction_gain_factor_rectification = ...
    interaction_gain_factor_rectification;
loss_spectrum_parameters.interaction_gain_factor_photodember =...
    interaction_gain_factor_photodember;
loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
    v_struct.(strcat('angle_',num2str(theta_pol_degree))) + ...
    interaction_gain_factor_photodember * circshift(v_struct.(strcat('photodember')),[delay 0]);

[~ , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);

            
end


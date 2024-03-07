function [eels] = setup_parameters_eels_photodember(pump_power_nj, laser_spot_size_fwhm)
%DEFAULT_PARAMETERS Summary of this function goes here
%   Detailed explanation goes here

% Laser parameters
% laser_parameters.pulse_energy_experiment = 10 * 1e-9;
laser_parameters.pulse_energy_experiment = 1 * pump_power_nj * 1e-9;
% laser_parameters.pulse_energy_gain_factor = 0.072;
laser_parameters.pulse_energy_gain_factor = 0.014;
laser_parameters.laser_spot_fwhm = laser_spot_size_fwhm;
laser_parameters.theta_pol = 90*pi/180;
laser_parameters.laser_pulse_time_fwhm = 50e-15;

% Discretization parameters
discretization_params.x0 = 0;
discretization_params.y0 = -1e-6;
discretization_params.ddt = 10e-15;   
% discretization_params.delay_max = 2.5e-12;
discretization_params.delay_max = 2 * 1.5e-12;
discretization_params.fs = 2.4e15;    
discretization_params.l = 2.4e4;
% discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
discretization_params.t0 = -0.5e-12;
% discretization_params.t0 = 0e-12;
discretization_params.xprime_max = 3 * Laser.calculate_sigma(laser_parameters.laser_spot_fwhm) ;
discretization_params.d_xprime = 4e-2 * 3 * Laser.calculate_sigma(laser_parameters.laser_spot_fwhm) ;
discretization_params.yprime_max = 1e-6;
discretization_params.d_yprime = 4e-2 * 1e-6;
discretization_params.zprime_max = discretization_params.xprime_max;
discretization_params.d_zprime = discretization_params.d_xprime;
discretization_params.ddz = 2e-6;     
discretization_params.zmax = 1e-4;
discretization_params.z_max = 30e-6;
% discretization_params.z_max = 90e-6;

% UTEM parameters
utem_parameters.electron_total_energy = 1.1;
utem_parameters.electron_total_energy = 0.94;
utem_parameters.electron_total_time_fs = 360;
utem_parameters.electron_time_coherent_fwhm_fs = 50;
utem_parameters.electron_theta = -7*pi/180;
utem_parameters.electron_velocity_c = 0.7;

% Numerical parameters
numerical_parameters.tc_subsampling = 30;
numerical_parameters.subsampling_factor = 60;

    
laser = Laser(laser_parameters);

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

% [w, e_w, t_w_store] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
%     discretization.energy, discretization.deltat);


eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
eels_parameters.materialmaterial.phase = 0;
eels = EELS(eels_parameters);

% tt = discretization.t;
% zz = discretization.z;

% params_photodember.eels = eels;


end

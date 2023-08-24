function [laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2(spot_size_fwhm)
%DEFAULT_PARAMETERS Summary of this function goes here
%   Detailed explanation goes here
laser_parameters.pulse_energy_experiment = 10 * 1e-9;
laser_parameters.pulse_energy_gain_factor = 0.014;
if nargin < 1
    laser_parameters.laser_spot_fwhm = 40e-6;
else
    laser_parameters.laser_spot_fwhm = spot_size_fwhm;
end

laser_parameters.theta_pol = 45*pi/180;
laser_parameters.laser_pulse_time_fwhm = 50e-15;

discretization_params.x0 = 0;
discretization_params.y0 = -1e-6;
discretization_params.ddt = 10e-15;   discretization_params.delay_max = 2.5e-12;
discretization_params.fs = 2.4e15;    discretization_params.l = 2.4e4;
discretization_params.t0 = -0.5e-12;


discretization_params.xprime_max = 3 * Laser.calculate_sigma(laser_parameters.laser_spot_fwhm) ;
discretization_params.d_xprime = 4e-2 * 3 * Laser.calculate_sigma(laser_parameters.laser_spot_fwhm) ;
discretization_params.d_xprime = 2e-6 ;


discretization_params.yprime_max = 1e-6;
discretization_params.d_yprime = 4e-2 * 1e-6;
discretization_params.zprime_max = discretization_params.xprime_max;
discretization_params.d_zprime = discretization_params.d_xprime;
discretization_params.ddz = 2e-6;     discretization_params.zmax = 1e-4;
discretization_params.z_max = 30e-6;

utem_parameters.electron_total_energy = 1.1;
utem_parameters.electron_total_time_fs = 360;
utem_parameters.electron_time_coherent_fwhm_fs = 50;
utem_parameters.electron_theta = -7*pi/180;
utem_parameters.electron_velocity_c = 0.7;

% % subsampling
numerical_parameters.tc_subsampling = 30;
numerical_parameters.subsampling_factor = 60;

end
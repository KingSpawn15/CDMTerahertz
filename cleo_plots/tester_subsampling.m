clearvars;

laser_parameters.pulse_energy_experiment = 1 * 1e-9;
laser_parameters.pulse_energy_gain_factor = 0.014;
laser_parameters.laser_spot_fwhm = 80e-6;
laser_parameters.theta_pol = 90*pi/180;
laser_parameters.laser_pulse_time_fwhm = 50e-15;

discretization_params.x0 = 0;
discretization_params.y0 = -1e-6;
discretization_params.ddt = 10e-15;   discretization_params.delay_max = 2.5e-12;
discretization_params.fs = 2.4e15;    discretization_params.l = 2.4e4;
discretization_params.t0 = -0.5e-12;
discretization_params.xprime_max = 3 * Laser.calculate_sigma(laser_parameters.laser_spot_fwhm);
discretization_params.d_xprime = 4e-2 * 3 * Laser.calculate_sigma(laser_parameters.laser_spot_fwhm);
discretization_params.yprime_max = 1e-6;
discretization_params.d_yprime = 4e-2 * 1e-6;
discretization_params.zprime_max = discretization_params.xprime_max;
discretization_params.d_zprime = discretization_params.d_xprime;
discretization_params.ddz = 2e-6;     discretization_params.zmax = 1e-4;
discretization_params.z_max = 60e-6;

utem_parameters.electron_total_energy = 0.94;
utem_parameters.electron_total_time_fs = 360;
utem_parameters.electron_time_coherent_fwhm_fs = 50;
utem_parameters.electron_theta = -7*pi/180;
utem_parameters.electron_velocity_c = 0.7;

% % subsampling
numerical_parameters.tc_subsampling = 15;
numerical_parameters.subsampling_factor = 60;



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
interaction_v_pd =  ChargeDynamics.interaction_potential_photodember(discretization, material,...
    laser , numerical_parameters);

% nteraction_potential_rectification(discretization, material,...
%                 laser , electron, numerical_parameters)
toc;

x
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
% method = "photodember";
% interaction_gain_factor_photodember = 31.2;
% loss_spectrum_parameters.method = method;
% 
% loss_spectrum_parameters.interaction_gain_factor_rectification = 1;
% loss_spectrum_parameters.interaction_gain_factor_photodember =...
%     interaction_gain_factor_photodember;
interaction_gain_factor_photodember = 0.07;
loss_spectrum_parameters.interact_v = interaction_v * interaction_gain_factor_photodember;

[psi_sub_com , psi_incoherent_com] = eels.energy_loss_spectrum(loss_spectrum_parameters);

% 
[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
        discretization.energy, discretization.deltat);
imagesc(e_w,t_w, psi_sub_com);
ylim([-1 , 1.5]);

colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.XTick = -4:2:4;
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
ax.YTick = -1:0.5:1.5;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
clearvars;


[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();


laser_parameters.pulse_energy_experiment = 1e-9;
discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
discretization_params.delay_max = 2 * 1.5e-12;

utem_parameters.electron_total_energy = 0.94;
laser_parameters.laser_pulse_time_fwhm = 650e-15;
laser_parameters.theta_pol = 90*pi/180;

laser = Laser(laser_parameters);

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

[w, e_w, t_w_store] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);


eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
eels_parameters.materialmaterial.phase = 0;
eels = EELS(eels_parameters);

loss_spectrum_parameters.method = 'rectification';
loss_spectrum_parameters.interaction_gain_factor_rectification = 1;
loss_spectrum_parameters.interaction_gain_factor_photodember = 0;
%%

interact_v_or_store = eels.interaction_v(loss_spectrum_parameters);

%%
loss_spectrum_parameters.method = 'photodember';
loss_spectrum_parameters.interaction_gain_factor_rectification = 0;
loss_spectrum_parameters.interaction_gain_factor_photodember = 1;


interact_v_pd_store = eels.interaction_v(loss_spectrum_parameters);


%%
alpha_pd_0 =  .07;
alpha_or_0 = 20 * 1.3 * 1.2;

interact_v_pd = circshift(interact_v_pd_store, [18,0]);
interact_v_or = circshift(interact_v_or_store, [-15 ,0]);


t_w = t_w_store-0.2;

alpha_pd = alpha_pd_0; alpha_or =  alpha_or_0;
loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
    interact_v_or * alpha_or;


% loss_spectrum_parameters.interact_v = interact_v_pd * 0.07*exp(1i*0);

[psi_sub_comb , psi_incoherent_comb] = eels.energy_loss_spectrum(loss_spectrum_parameters);


alpha_pd = 0; alpha_or =  alpha_or_0;
loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
    interact_v_or * alpha_or;


% loss_spectrum_parameters.interact_v = interact_v_pd * 0.07*exp(1i*0);

[psi_sub_or , psi_incoherent_or] = eels.energy_loss_spectrum(loss_spectrum_parameters);

alpha_pd = alpha_pd_0; alpha_or =  0;
loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
    interact_v_or * alpha_or;


% loss_spectrum_parameters.interact_v = interact_v_pd * 0.07*exp(1i*0);

[psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);

imagesc(e_w,t_w, psi_incoherent_comb);

ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
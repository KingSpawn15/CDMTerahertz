clearvars;


[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();

laser = Laser(laser_parameters);

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);


eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
eels = EELS(eels_parameters);

loss_spectrum_parameters.method = 'photodember';
loss_spectrum_parameters.interaction_gain_factor_rectification = 0;
loss_spectrum_parameters.interaction_gain_factor_photodember = 1;
interact_v = eels.interaction_v(loss_spectrum_parameters);


loss_spectrum_parameters.interact_v = interact_v * 0.07*exp(1i*0);

[psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);

imagesc(e_w,t_w, psi_incoherent);

ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);



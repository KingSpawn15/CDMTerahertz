clearvars;

%%
for pangle = [0, 45, 135, 90]
export_dir = 'examples/phase_study/original/revision/';
[status, msg, msgID] = mkdir(export_dir);
base_filename = 'eels_pd_electron_';
pol_angle = pangle;
export_file_name = strcat(export_dir,base_filename,num2str(pol_angle),'.png');

%%
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();
laser = Laser(laser_parameters);

discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;
discretization_params.delay_max = 1.5e-12;

utem_parameters.electron_total_energy = 0.94;

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);



eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
% laser_parameters.laser_pulse_time_fwhm = 500e-15;
laser_parameters.theta_pol = pol_angle * pi / 180;
laser = Laser(laser_parameters);
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
% eels_parameters.material.phase = 0;
% eels_parameters.material.gamma_factor = 0.8;
eels = EELS(eels_parameters);

loss_spectrum_parameters.method = 'photodember';
loss_spectrum_parameters.interaction_gain_factor_rectification = 0;
loss_spectrum_parameters.interaction_gain_factor_photodember = 1;
interact_v_pd = eels.interaction_v(loss_spectrum_parameters);

loss_spectrum_parameters.method = 'rectification';
loss_spectrum_parameters.interaction_gain_factor_rectification = 1;
loss_spectrum_parameters.interaction_gain_factor_photodember = 0;
interact_v_or = eels.interaction_v(loss_spectrum_parameters);

interact_v_pd_store = interact_v_pd;
interact_v_or_store = interact_v_or;
%%

% alpha_pd_0 =  0.001;
% alpha_or_0 = 4*exp(1i*pi/180*0);

alpha_pd_0 =  0.07;
alpha_or_0 = 0.2;
% alpha_or_0 = 0;

interact_v_pd = circshift(interact_v_pd_store, [-10,0]);
interact_v_or = circshift(interact_v_or_store, [-15,0]);
alpha_pd = alpha_pd_0; alpha_or =  0;
loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
    interact_v_or * alpha_or;
[psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);

alpha_pd = 0; alpha_or =  alpha_or_0;
loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
    interact_v_or * alpha_or;
[psi_sub_or , psi_incoherent_or] = eels.energy_loss_spectrum(loss_spectrum_parameters);

alpha_pd = alpha_pd_0; alpha_or =  alpha_or_0;
loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
    interact_v_or * alpha_or;
[psi_sub_com , psi_incoherent_com] = eels.energy_loss_spectrum(loss_spectrum_parameters);

close all;
tiledlayout('flow');
nexttile;
imagesc(e_w,t_w, psi_sub_pd);
ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ax.XTick = -4:2:4;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

nexttile;
imagesc(e_w,t_w, psi_sub_or);
ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ax.XTick = -4:2:4;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

nexttile;
imagesc(e_w,t_w, psi_sub_com);
ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ax.XTick = -4:2:4;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

nexttile;
imagesc(e_w,t_w, psi_incoherent_pd);
ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ax.XTick = -4:2:4;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

nexttile;
imagesc(e_w,t_w, psi_incoherent_or);
ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.XTick = -4:2:4;
ax.YTick = -1:0.5:1.5;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

nexttile;
imagesc(e_w,t_w, psi_incoherent_com);
ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ax.XTick = -4:2:4;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

set(gcf,'Position',[100, 50, 1350, 450*2]);
exportgraphics(gcf, export_file_name,'resolution' , 400);
end
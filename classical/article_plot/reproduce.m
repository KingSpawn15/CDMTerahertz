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

imagesc(e_w,t_w, psi_sub_pd);

ylim([-1 , 1.5]);
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:1.5;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

%%

zz = -1e-6;
dd = .5e-3;
lambda = 800e-9;
tau = 30e-15;
sigma_z = 80e-6;
[t0_vec, eels_calc] = eels_theoretical(tau, lambda, dd, zz, sigma_z);
eels_t = interp1(t0_vec.',eels_calc.',t_w,'linear','extrap');



[~, eels_ind_pd] = max(psi_sub_pd,[],2);
[~, eels_ind_or] = max(psi_sub_or,[],2);
eels_ind_pd = circshift(eels_ind_pd,0);
e_exc_pd = e_w(eels_ind_pd);
e_exc_or = e_w(eels_ind_or);

e_exc_pd = circshift(e_exc_pd,15)

factor_th = 0.9;
psi_assemb_pd  = assemble_psi_sub(t_w, e_w, e_exc_pd, psi_sub_or);
psi_assemb_or  = assemble_psi_sub(t_w, e_w, e_exc_or, psi_sub_or);
psi_assemb_comb  = assemble_psi_sub(t_w, e_w, e_exc_or + e_exc_pd, psi_sub_or);
psi_assemb_theory  = assemble_psi_sub(t_w, e_w, eels_t * factor_th, psi_sub_or);
psi_assemb_theory_comb  = assemble_psi_sub(t_w, e_w, e_exc_pd + eels_t * factor_th, psi_sub_or);
psi_assemb_theory_neg  = assemble_psi_sub(t_w, e_w, e_exc_pd - eels_t * factor_th, psi_sub_or);
close all
plo1 = plot(t_w, factor_th * eels_t, t_w, e_exc_or);
plo1(1).LineWidth = 2;
plo1(2).LineWidth = 1;
xlim([-1,1.5]);
% xlim([-1,1.5])
exportgraphics(gcf,'classical/article_plot/comparison.png','Resolution',500)

%%

% psi_incoherent_assemb = eels.incoherent_convolution(psi_assemb_comb , w, t_w, e_w);
% imagesc(e_w,t_w, psi_incoherent_assemb);
% ylim([-1 , 1.5]);
close all

figure;
psi_incoherent_comb = eels.incoherent_convolution(psi_assemb_theory_comb, w, t_w, e_w);
imagesc(e_w,t_w, psi_incoherent_comb);
ylim([-1 , 1.5]);



colormap jet
axis square
ax = gca;
ax.FontSize = 20;
ax.FontName = 'helvetica'
ax.LineWidth = 1;
ax.YTick = -.5:0.5:1.5;
ax.XTick = -4:2:4;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',22, 'FontName' , 'helvetica');
xlabel('Energy [eV]','Color',[0.3 0.3 0.3],'FontSize',22, 'FontName' , 'helvetica');
exportgraphics(gcf,'classical/article_plot/comparison_90.png','Resolution',500)


figure;
psi_incoherent_neg = eels.incoherent_convolution(psi_assemb_theory_neg, w, t_w, e_w);
imagesc(e_w,t_w, psi_incoherent_neg);
ylim([-1 , 1.5]);

colormap jet
axis square
ax = gca;
ax.FontSize = 20;
ax.FontName = 'helvetica'
ax.LineWidth = 1;
ax.YTick = -.5:0.5:1.5;
ax.XTick = -4:2:4;
ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',22, 'FontName' , 'helvetica');
xlabel('Energy [eV]','Color',[0.3 0.3 0.3],'FontSize',22, 'FontName' , 'helvetica');
exportgraphics(gcf,'classical/article_plot/comparison_0.png','Resolution',500)
%%
function psi_assemb = assemble_psi_sub(t_w, e_w, eels, psi_sub)
    psi_assemb = 1;

    eels_w = eels.';
    eels_w = repmat(eels_w,[1,length(e_w)]);

    e_w_mat = e_w;
    e_w_mat = abs(repmat(e_w_mat,[length(t_w),1]) - eels_w);
    
    [~, b] = min(e_w_mat,[],2);
    
    psi_assemb = psi_sub * 0;
    row = 1;
    for col = b.'
        psi_assemb(row,col) = 1;
        row = row + 1;
    end

    

end
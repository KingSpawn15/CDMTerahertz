clearvars -except TPD ZPD EPD psi_sub_pd  psi_incoherent_pd
close all
% delete any existing parallel pool
% delete(gcp('nocreate'))
% % create a new parallel pool with 6 workers
% parpool(6);
% parfor ii = 1:1024
% end

setdir = 'investigation/noshift/';
or_spot_sigma = 50e-6;
pd_spot_fwhm = 40e-6;
pd_z_max = 90e-6;
pulse_energy_experiment_nj = 10;



offset = 101;
theta = @(offset,angle) mod(fix(2*(-offset + angle)),360) - 180;
angle_0 = 56 - 45;
angle_45 = 56 - 22;
angle_90 = 56;
    
[com_0, ~, ~, eels_measure_0] = get_measurement_3(angle_0);
[com_90, ~, ~, eels_measure_90] = get_measurement_3(angle_90);
[com_45, deltat, energy, eels_measure_45] = get_measurement_3(angle_45);
com_or_90_hypo = com_90 - com_45;

% pulse_energy = 10;
[eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max);
% interact_v_pd_store = utils.eels_pdor(eels,'photodember');
% interact_v_pd = circshift(interact_v_pd_store, [18,0]);


% alpha_pd_0 =  .07;
% loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0;
% [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);
% [~, eels_ind_pd] = max(psi_sub_pd,[],2);
% eels_ind_pd = circshift(eels_ind_pd,0);
% e_exc_pd = e_w(eels_ind_pd);


[TPD, ZPD, EPD,psi_sub_pd , psi_incoherent_pd] = electric_field_photodember(pulse_energy_experiment_nj, ...
    pd_spot_fwhm, pd_z_max);
[TOR, ZOR, EOR] = electric_field_rectification(or_spot_sigma);
EPD = EPD*(1.26/6.34);
%%

%%

EOR1 = EOR(ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:)),TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:)));
TOR1 = TOR(TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:)),ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:)));
ZOR1 = ZOR(TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:)),ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:)));
EORintrap = interp2(TOR1.', ZOR1.', EOR1, TPD.', ZPD.', 'linear', 0);
EPDintrap = interp2(TPD.', ZPD.', EPD, TOR.', ZOR.', 'linear', 0);

close all

[~,~] = mkdir(setdir);
FontSize = 10;
FontName = 'ariel';

% normalize = @(EE) real(EE)/ max(abs(real(EE(:))));
%%

[~, eels_ind_pd] = max(psi_sub_pd,[],2);
e_exc_pd = e_w(eels_ind_pd);
velec = 0.7 * 3 * 10^(8-12);
[~, ~, ~, t0_vec_pd, eels_cc_pd] = eels_pd(EPD.', tt * 1e12, zz , velec);
%%
tic;
[t0_vec_or, eels_cc_or] = eels_calc_or_test(EOR.', TOR, ZOR* 1e-6);
toc;
%%
tic;
[t0_vec_or, eels_cc_or] = utils_spectrum.calculate_spectrum_from_fields(EOR.', TOR, ZOR* 1e-6);
toc;
%%
[t0_vec_pd2, eels_cc_pd2] = utils_spectrum.calculate_spectrum_from_fields(EPDintrap.' , TOR, ZOR* 1e-6);
[t0_vec_pd2, eels_cc_comb_0] = utils_spectrum.calculate_spectrum_from_fields((EPDintrap - EOR).' , TOR, ZOR* 1e-6);
[t0_vec_pd2, eels_cc_comb_90] = utils_spectrum.calculate_spectrum_from_fields((EPDintrap + EOR).' , TOR, ZOR* 1e-6);
%%
close all
FontName = 'ariel';
FontSize = 10;
figure;
plot_tile(energy, deltat, eels_measure_90);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 150]);
hold on;
% plot(e_exc_pd , t_w);
plot(eels_cc_or + eels_cc_pd2, t0_vec_pd2);
set(gca,'ydir','reverse')

figure;
plot_tile(e_w,t_w, psi_incoherent_comb_90);
hold on;
plot(eels_cc_or + eels_cc_pd2, t0_vec_pd2 );
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 150]);
%%
figure;
imagesc(TPD(:,1), ZPD(1,:), EPD);
set(gca,'FontSize',FontSize);
xlim([-.3,1.5]);
ylim([-100,100]);
xticks(-.3:.3:1.5)
colormap(utils.redblue);
pbaspect([2 1 1])
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
colorbar;
exportgraphics(gcf, [setdir, 'field_photodember.png'],'resolution', 300);

figure;
imagesc(TOR1(:,1), ZOR1(1,:), EOR1);
set(gca,'FontSize',FontSize);
xlim([-.3,1.5]);
xlim([0,1.5]);
ylim([-100,100]);
xticks(-.3:.3:1.5)
colormap(utils.redblue);
pbaspect([2 1 1])
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
colorbar;
exportgraphics(gcf, [setdir, 'field_rectification.png'],'resolution', 300);

figure;
imagesc(TPD(:,1), ZPD(1,:), EORintrap + EPD);
set(gca,'FontSize',FontSize);
xlim([-.3,1.5]);
xlim([0,1.5]);
ylim([-100,100]);
xticks(-.3:.3:1.5)
colormap(utils.redblue);
pbaspect([2 1 1])
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
colorbar;
exportgraphics(gcf, [setdir, 'field_combined.png'],'resolution', 300);

%%
image_name = 'power_mod';
sigma_z = or_spot_sigma;
factor_th = 1;
% shift_tw = 0.15;
shift_tw = 0;
t_w = t_w_store-shift_tw; 

for shift = [0]

    psi_or = rectification_eels(t_w, sigma_z, shift);
    psi_incoherent_comb_90 = combined_eels(psi_sub_pd, factor_th * psi_or, eels, w, t_w, e_w);
    psi_incoherent_comb_0 = combined_eels(psi_sub_pd, -factor_th * psi_or, eels, w, t_w, e_w);
    psi_incoherent_or_90 = calculate_eels_or(factor_th * psi_or, eels, w, t_w, e_w);
    psi_incoherent_or_0 = calculate_eels_or(-factor_th * psi_or, eels, w, t_w, e_w);
    
    
    figure;
    FontName = 'ariel';
    FontSize = 10;
    ttt = tiledlayout(2,4,"TileSpacing","compact");
    ttt.Padding = "loose";
    
    plot_tile(energy, deltat, eels_measure_0);
    set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3]);
    plot_tile(e_w,t_w, psi_incoherent_comb_0);
    plot_tile(e_w,t_w, psi_incoherent_pd);
    plot_tile(e_w,t_w, psi_incoherent_or_0);
    plot_tile(energy, deltat, eels_measure_90);
    set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
    plot_tile(e_w,t_w, psi_incoherent_comb_90);
    set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
    plot_tile(e_w,t_w, psi_incoherent_pd);
    plot_tile(e_w,t_w, psi_incoherent_or_90);

set(gcf,'Position',[200,200,200 + 2 * 300,200 +  300]); 


exportgraphics(gcf, [setdir,image_name,num2str(shift),'.png'], 'Resolution',300);


end
    
%%


% function [t0_vec, eels] = eels_calc_or_test(EOR, TOR, ZOR)
% sigma_z = 50e-6;
% % vec_z = (-20 : 0.05 : 20 ).'*sigma_z;
% % [T, Z] = ndgrid(t, vec_z);
% % ET = repmat(ethz, 1, length(vec_z)) .* exp(-Z.^2 / (2 * sigma_z^2));
% % ET(abs(Z)>100e-6) =0;
% 
% velec = 0.7 * 3 * 10^(8-12);
% 
% 
% 
% t0_vec = (-5:0.02:5);
% eels = zeros(length(t0_vec),1);
% ind = 1;
% for t0 = t0_vec
% 
%     zz_l = (-5 : 0.1: 5).' .* sigma_z;
%     tt_l = zz_l / velec + t0;
%     et_l = interp2(TOR', ZOR', EOR', tt_l, zz_l,"linear",0);
% %     et_l(isnan(et_l)) = 0;
%     eels(ind) = trapz(zz_l,et_l);
%     ind = ind + 1;
% end
% 
% 
% t0_vec = t0_vec.';
% 
% 
% 
% 
% end

% function [T, Z, ET, t0_vec, eels] = eels_pd(ET, t, z , velec)
%     
%     zmax = 2.5e-4;
%     [T, Z] = ndgrid(t, z);
% 
%     t0_vec = (-2:0.02:2);
%     eels = zeros(length(t0_vec),1);
%     ind = 1;
%     for t0 = t0_vec
%     
%         zz_l = (-zmax : 0.1 * zmax: zmax).';
%         tt_l = zz_l / velec + t0;
%         et_l = interp2(T', Z', ET', tt_l, zz_l,"linear",0);
%     %     et_l(isnan(et_l)) = 0;
%         eels(ind) = trapz(zz_l,et_l);
%         ind = ind + 1;
%     end
%     
%     
%     t0_vec = t0_vec.';
% 
% 
% end

function [T, Z, ET] = electric_field_rectification(sigma_z)

    zz = -1e-6;
    dd = .5e-3;
    lambda = 800e-9;
    tau = 30e-15;
%     sigma_z = 50e-6;
    [~, ~,T, Z, ET] = eels_theoretical_2(tau, lambda, dd, zz, sigma_z);
    ET = ET.';
    Z = Z * 1e6;

end


function [TT, ZZ, electric_field_zt,psi_sub_pd , psi_incoherent_pd] = electric_field_photodember(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max)
% This function returns tt, zz, and electric_field_zt

%     % delete any existing parallel pool
%     delete(gcp('nocreate'))
%     
%     % create a new parallel pool with 6 workers
%     parpool(6);
%     
    % set up the eels parameters
    [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max);
    
    % calculate the photodember interaction voltage
    interact_v_pd_store = utils.eels_pdor(eels,'photodember');
    
    % shift the interaction voltage by 18 steps along the first dimension
    interact_v_pd = circshift(interact_v_pd_store, [0,0]);
    
    % adjust the time window by subtracting 0.2
%     t_w = t_w_store-0.2;
    
    % set the photodember coefficient to 0.07
    alpha_pd_0 =  0.05;
    
    % multiply the interaction voltage by the coefficient
    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 ;
    
    % calculate the energy loss spectrum for photodember
    [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);
    
    % find the maximum of the energy loss spectrum along the second dimension
%     [~, eels_ind_pd] = max(psi_sub_pd,[],2);
    
    % find the corresponding excitation energy
%     e_exc_pd = e_w(eels_ind_pd);
    
    % create a 2d grid of tt and zz
    [TT, ZZ] = ndgrid(tt * 1e12,zz * 1e6);
    
    % use the derivative_f_dz function to calculate the electric field along z and t
    electric_field_zt = derivative_f_dz(interact_v_pd_store, zz, tt);
    
    % replace any NaN values with zero
    electric_field_zt(isnan(electric_field_zt)) = 0;

end

function plot_tile(x, y, z)
    FontSize = 10;
    FontName = 'ariel';
    nexttile
    imagesc(x, y, z);
    ylim([-1,1.5]);
    xlim([-5,5]);
    colormap jet
    axis square
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3]);
end

function [com, deltat, energy, eels_measure] = get_measurement_3(angle)


    [eels_measure, energy,...
        time] = measurement_plots.data_measurement(angle);
    
    deltat = time(time > -1 & time < 1.5);
    eels_measure = eels_measure(time > -1 & time < 1.5,:);
    
    x_c = zeros(1, size(deltat,2));
    
    eels_measure = eels_measure ./ sqrt(sum(eels_measure.^2, 2));
    
    
    for i = 1:size(eels_measure,1)
        [~, ind_max] = max(eels_measure(i,:));
        x_c(i) = ind_max;
    end
    
    com = energy(floor(x_c));
end


function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end

function [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max)
%     spot_fwhm = 40e-6;
%     [laser_parameters,discretization_params, utem_parameters,...
%         numerical_parameters] = default_parameters_2(spot_fwhm);

    
%     if nargin < 1
%         pulse_energy_experiment_nj = 10;
%     end

     [laser_parameters,discretization_params, utem_parameters,...
        numerical_parameters] = default_parameters_2(pd_spot_fwhm);

    laser_parameters.pulse_energy_experiment = 0.1 * pulse_energy_experiment_nj *1e-9;
    discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
    discretization_params.delay_max = 2 * 1.5e-12;
    discretization_params.z_max = pd_z_max;
    
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

    tt = discretization.t;
    zz = discretization.z;
end

function psi_assemb = assemble_psi_sub(t_w, e_w, eels)

    eels_w = eels.';
    eels_w = repmat(eels_w,[1,length(e_w)]);

    e_w_mat = e_w;
    e_w_mat = abs(repmat(e_w_mat,[length(t_w),1]) - eels_w);
    
    [~, b] = min(e_w_mat,[],2);
    
    psi_assemb = zeros(length(t_w), length(e_w));
    row = 1;
    for col = b.'
        psi_assemb(row,col) = 1;
        row = row + 1;
    end

    

end


function psi_or = rectification_eels(t_w, sigma_z, shift)
      
    zz = -1e-6;
    dd = .5e-3;
    lambda = 800e-9;
    tau = 30e-15;
    [t0_vec, eels_calc_cc,~, ~, ~] = eels_theoretical_2(tau, lambda, dd, zz, sigma_z);
    eels_calc_cc = circshift(eels_calc_cc, shift); 
    psi_or = interp1(t0_vec.',eels_calc_cc.',t_w,'linear','extrap');

end

function psi_incoherent_comb = combined_eels(psi_pd, psi_or, eels_obj, w, t_w, e_w)
    
    [~, eels_ind_pd] = max(psi_pd,[],2);
    psi_pd_extracted = e_w(eels_ind_pd);

    psi = psi_or + psi_pd_extracted;
    psi_comb  = assemble_psi_sub(t_w, e_w, psi);
    psi_incoherent_comb = eels_obj.incoherent_convolution(psi_comb, w, t_w, e_w);

end

function psi_incoherent_comb = calculate_eels_or(psi_or, eels_obj, w, t_w, e_w)
    
   
    psi = psi_or ;
    psi_comb  = assemble_psi_sub(t_w, e_w, psi);
    psi_incoherent_comb = eels_obj.incoherent_convolution(psi_comb, w, t_w, e_w);

end

function dfdz = derivative_f_dz(f_values, z_values, t_values)
    % Check if the input matrices have compatible sizes
    if size(f_values, 1) ~= length(z_values) || size(f_values, 2) ~= length(t_values)
        error('Input matrix dimensions do not match with vector sizes.');
    end

    % Calculate the derivative using the finite difference method
    dz = z_values(2) - z_values(1); % Assuming uniform grid spacing
    dfdz = diff(f_values, 1, 1) ./ dz;

    % Pad the dfdz matrix with NaNs to maintain the same size as f_values
    dfdz = [dfdz; NaN(1, size(dfdz, 2))];
end


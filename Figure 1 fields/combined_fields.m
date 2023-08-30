clear all
close all

[TPD, ZPD, EPD] = electric_field_photodember();
[TOR, ZOR, EOR] = electric_field_rectification();

EOR1 = EOR(ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:)),TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:)));
TOR1 = TOR(TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:)),ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:)));
ZOR1 = ZOR(TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:)),ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:)));
EORintrap = interp2(TOR1.', ZOR1.', EOR1, TPD.', ZPD.', 'linear', 0);


setdir = 'Figure 1 fields/results/';
FontSize = 10;
FontName = 'ariel';

normalize = @(EE) real(EE)/ max(abs(real(EE(:))));

%%
figure;
imagesc(TPD(:,1), ZPD(1,:), normalize(EPD),[-1,1]);
set(gca,'FontSize',FontSize);
xlim([-.3,1.5]);
ylim([-100,100]);
xticks(-.3:.3:1.5)
colormap(utils.redblue);
pbaspect([2 1 1])
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
exportgraphics(gcf, [setdir, 'field_photodember.png'],'resolution', 300);

figure;
imagesc(TOR1(:,1), ZOR1(1,:), normalize(EOR1),[-1,1]);
set(gca,'FontSize',FontSize);
xlim([-.3,1.5]);
ylim([-100,100]);
xticks(-.3:.3:1.5)
colormap(utils.redblue);
pbaspect([2 1 1])
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
exportgraphics(gcf, [setdir, 'field_rectification.png'],'resolution', 300);

figure;
imagesc(TPD(:,1), ZPD(1,:), normalize(EORintrap) + normalize(EPD),[-1,1]);
set(gca,'FontSize',FontSize);
xlim([-.3,1.5]);
ylim([-100,100]);
xticks(-.3:.3:1.5)
colormap(utils.redblue);
pbaspect([2 1 1])
set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
exportgraphics(gcf, [setdir, 'field_combined.png'],'resolution', 300);


%% Function to set axis properties

function [T, Z, ET] = electric_field_rectification()

    zz = -1e-6;
    dd = .5e-3;
    lambda = 800e-9;
    tau = 30e-15;
    sigma_z = 50e-6;
    [~, ~,T, Z, ET] = eels_theoretical_2(tau, lambda, dd, zz, sigma_z);
    ET = ET.';
    Z = Z * 1e6;

end


function [TT, ZZ, electric_field_zt] = electric_field_photodember()
% This function returns tt, zz, and electric_field_zt

    % delete any existing parallel pool
    delete(gcp('nocreate'))
    
    % create a new parallel pool with 6 workers
    parpool(6);
    
    % set up the eels parameters
    [eels, w, e_w, t_w_store, tt, zz]= eels_setup();
    
    % calculate the photodember interaction voltage
    interact_v_pd_store = utils.eels_pdor(eels,'photodember');
    
    % shift the interaction voltage by 18 steps along the first dimension
    interact_v_pd = circshift(interact_v_pd_store, [18,0]);
    
    % adjust the time window by subtracting 0.2
    t_w = t_w_store-0.2;
    
    % set the photodember coefficient to 0.07
    alpha_pd_0 =  0.07;
    
    % multiply the interaction voltage by the coefficient
    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 ;
    
    % calculate the energy loss spectrum for photodember
    [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);
    
    % find the maximum of the energy loss spectrum along the second dimension
    [~, eels_ind_pd] = max(psi_sub_pd,[],2);
    
    % find the corresponding excitation energy
    e_exc_pd = e_w(eels_ind_pd);
    
    % create a 2d grid of tt and zz
    [TT, ZZ] = ndgrid(tt * 1e12,zz * 1e6);
    
    % use the derivative_f_dz function to calculate the electric field along z and t
    electric_field_zt = derivative_f_dz(interact_v_pd_store, zz, tt);
    
    % replace any NaN values with zero
    electric_field_zt(isnan(electric_field_zt)) = 0;

end

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,YTickLabel,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
    ax.LineWidth = LineWidth;
    ax.YTick = YTick;
    if ~isempty(YTickLabel)
        ax.XTick = YTickLabel;
    end
    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end


function [T, Z, ET, t0_vec, eels] = eels_pd(ET, t, z , velec)
    
    zmax = 0.5e-4;
    [T, Z] = ndgrid(t, z);

    t0_vec = (-2:0.02:2);
    eels = zeros(length(t0_vec),1);
    ind = 1;
    for t0 = t0_vec
    
        zz_l = (-zmax : 0.1 * zmax: zmax).';
        tt_l = zz_l / velec + t0;
        et_l = interp2(T', Z', ET', tt_l, zz_l);
    %     et_l(isnan(et_l)) = 0;
        eels(ind) = trapz(zz_l,et_l);
        ind = ind + 1;
    end
    
    
    t0_vec = t0_vec.';


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

function dfdt = derivative_f_dt(f, z, t)
% f is a 2d matrix of size m x n, where m is the length of z and n is the length of t
% z is a vector of size m x 1, containing the values of z
% t is a vector of size n x 1, containing the values of t
% ft is a 2d matrix of size m x n, containing the values of ft/dt(z,t)

% initialize ft with zeros
dfdt = zeros(size(f));

% calculate the constant spacing between t values
dt = t(2) - t(1);

% use central difference formula to approximate ft/dt(z,t) for each value of t, except the endpoints
dfdt(:,2:end-1) = (f(:,3:end) - f(:,1:end-2)) / (2 * dt);

% use forward difference formula to approximate ft/dt(z,t) at the first endpoint of t
dfdt(:,1) = (f(:,2) - f(:,1)) / dt;

% use backward difference formula to approximate ft/dt(z,t) at the last endpoint of t
dfdt(:,end) = (f(:,end) - f(:,end-1)) / dt;

end

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


function [eels, w, e_w, t_w_store, tt, zz]= eels_setup()
%     spot_fwhm = 40e-6;
%     [laser_parameters,discretization_params, utem_parameters,...
%         numerical_parameters] = default_parameters_2(spot_fwhm);
     [laser_parameters,discretization_params, utem_parameters,...
        numerical_parameters] = default_parameters_2();
    
    laser_parameters.pulse_energy_experiment = 1e-9;
    discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
    discretization_params.delay_max = 2 * 1.5e-12;
    discretization_params.z_max = 30e-6;
    
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

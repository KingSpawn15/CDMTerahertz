clearvars;


[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();


laser_parameters.pulse_energy_experiment = 1e-9;
discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
discretization_params.delay_max = 2 * 1.5e-12;
discretization_params.z_max = 100e-6;

utem_parameters.electron_total_energy = 0.94;
laser_parameters.laser_pulse_time_fwhm = 650e-15;
laser_parameters.theta_pol = 90*pi/180;
laser_parameters.laser_spot_fwhm = 40e-6;
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

% loss_spectrum_parameters.method = 'rectification';
% loss_spectrum_parameters.interaction_gain_factor_rectification = 1;
% loss_spectrum_parameters.interaction_gain_factor_photodember = 0;
% %%
% 
% interact_v_or_store = eels.interaction_v(loss_spectrum_parameters);


loss_spectrum_parameters.method = 'photodember';
loss_spectrum_parameters.interaction_gain_factor_rectification = 0;
loss_spectrum_parameters.interaction_gain_factor_photodember = 1;
%%
interact_v_pd_store = eels.interaction_v(loss_spectrum_parameters);


%%


alpha_pd_0 =  .07;
% alpha_or_0 = 20 * 1.3 * 1.2;

interact_v_pd = circshift(interact_v_pd_store, [18,0]);
% interact_v_or = circshift(interact_v_or_store, [-15 ,0]);


t_w = t_w_store-0.2;


alpha_pd = alpha_pd_0; 
loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd ;

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


close all
[~, eels_ind_pd] = max(psi_sub_pd,[],2);
e_exc_pd = e_w(eels_ind_pd);



% close all;
tt = discretization.t(:);
zz = discretization.z(:);

[TT, ZZ] = ndgrid(tt,zz);
vv = interact_v_pd_store;
electric_field_zt = derivative_f_dz(vv, zz, tt);
electric_field_zt(isnan(electric_field_zt)) = 0;

velec = 0.7 * 3 * 10^(8-12);

[T, Z, ET, t0_vec, eels_cc] = eels_pd(electric_field_zt.', tt * 1e12, zz , velec);
% figure;
plot(e_exc_pd, t_w ,'LineWidth',2)
ylim([-1,1.5]);
xlim([-5,5]);
hold on;
plot(eels_cc/max(eels_cc(:))*max(e_exc_pd(:)),t0_vec(:)  - 0.05,'r--','LineWidth',2)
legend('CDEM','Classical','FontSize',18)
hold off
set(gca,'YDir','reverse') ;
set(gca,'FontSize',14);
xticks(-4:2:4);
setdir = 'classical/results/';
exportgraphics(gcf,strcat(setdir,'photodember_classical_comaparison.png'),'resolution',400)

figure;
s=surf(ZZ, TT, real(electric_field_zt).' / max(real(electric_field_zt(:))));
s.EdgeColor = 'none';
%%

close all
figure;
imagesc(tt * 1e12, zz * 1e6, real(electric_field_zt)/ max(abs(real(electric_field_zt(:)))),[-1,1]);
set(gca,'FontSize',12);
% ylim([-120,120]);
xlim([-.3,1.5]);
xticks(-.3:.3:1.5)
colormap(redblue);
pbaspect([1 0.5 1])
colorbar;

exportgraphics(gcf,strcat(setdir,'photodember_e_phi_dx.png'),'resolution',400)


%%

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

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 
end


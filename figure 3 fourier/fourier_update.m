clear all
close all

delete(gcp('nocreate'))
% create a new parallel pool with 6 workers
parpool(6);

smoothing = @(matrix, window)  movmean(movmean(matrix, window, 1) , window, 2);

%% get measurements
angle = sort(0:2:180);
mean_eels_measure =  cell(1, length(angle));
eels_measure = cell(1, length(angle));
errs = cell(1, length(angle));
% [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
for ii = 1:length(angle)
[mean_eels_measure{ii}, deltat, energy, eels_measure{ii}, errs{ii}] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(angle(ii)));
end
mean_eels_exp = smoothing(cell2mat(mean_eels_measure.'),4);
errs_eels_exp = cell2mat(errs.');
%% build experiment spectra

spectra_matrix_exp = cell(1, length(angle));
for ii = 1:length(angle)
    [f_exp, spectra_matrix_exp{ii}] = builder_spectra(eels_measure{ii}, energy, deltat(:));
end
spectra_exp = smoothing(cell2mat(spectra_matrix_exp).',4);
spectra_exp = real(spectra_exp) ./ max(abs(real(spectra_exp(:))));

%% problem parameters
or_spot_sigma = 50e-6;
pd_spot_fwhm = 40e-6;
pd_z_max = 90e-6;
pulse_energy_experiment_nj = 10;

%% Photodember field
delete(gcp('nocreate'))
% create a new parallel pool with 6 workers
parpool(6);
[TPD, ZPD, EPD,psi_sub_pd, psi_incoherent_pd, ...
    eels, w, e_w, t_w, tt, zz] = electric_field_photodember(pulse_energy_experiment_nj, ...
    pd_spot_fwhm, pd_z_max);
EPD = EPD*(1.26/6.34);

%% Rectification field
rectification_param_z0 = -1e-6;
rectification_param_thickness = .5e-3;
rectification_param_wavelength = 800e-9;
rectification_param_pulse_time_tau = 30e-15;
rectification_param_pulse_laser_spot_sigma = or_spot_sigma;

[~, ~,TOR, ZOR, EOR] = eels_theoretical_2(rectification_param_pulse_time_tau, ...
    rectification_param_wavelength, ...
    rectification_param_thickness, ...
    rectification_param_z0, ...
    rectification_param_pulse_laser_spot_sigma);

EOR = EOR.';
ZOR = ZOR .* 1e6;

% Define the conditions
condition1 = ZOR(1,:) > min(ZPD(:)) & ZOR(1,:) < max(ZPD(:));
condition2 = TOR(:,1) > min(TPD(:)) & TOR(:,1) < max(TPD(:));

% Apply conditions to the arrays
EOR1 = EOR(condition1, condition2);
TOR1 = TOR(condition2, condition1);
ZOR1 = ZOR(condition2, condition1);
EORintrap = interp2(TOR1.', ZOR1.', EOR1, TPD.', ZPD.', 'linear', 0);
EPDintrap = interp2(TPD.', ZPD.', EPD, TOR.', ZOR.', 'linear', 0);
ZOR = ZOR .* 1e-6;
%% calculate base eels
kfactor = 1.2;
[~,eels_pd] = utils_spectrum.calculate_spectrum_from_fields(0 * EOR + ...
        kfactor * EPDintrap, TOR, ZOR);
[t0_vec,eels_or] = utils_spectrum.calculate_spectrum_from_fields(1 * EOR + ...
        0 * EPDintrap, TOR, ZOR);
%% combine eels to get max and spectra
mean_eels_theory = cell(1, length(angle));
psi_incoherent = cell(1, length(angle));
for ii = 1:length(angle)
    disp(angle(ii))
    mean_eels_theory{ii} = -cos(2 * angle(ii) * pi / 180) * eels_or + eels_pd;
    psi_incoherent{ii} =  eels.incoherent_convolution(utils_spectrum.spectrum_to_coherent_eels(t_w, e_w, mean_eels_theory{ii}, t0_vec),...
        w, t_w, e_w);
end
mean_eels_sim = smoothing(cell2mat(mean_eels_theory).',4);

%% build theory spectra
spectra_matrix_theory= cell(1, length(angle));
for ii = 1:length(angle)
    [f_theory, spectra_matrix_theory{ii}] = builder_spectra(psi_incoherent{ii}, e_w, t_w);
end
spectra_theory = smoothing(cell2mat(spectra_matrix_theory).',4);
spectra_theory = real(spectra_theory) ./ max(abs(real(spectra_theory(:))));

%% plot
close all
setdir = 'figure 3 fourier/results/';
FontName = 'ariel';
FontSize = 14;
dt_arr = -0.1;

% plot max eels experiment
figure;
image_name = 'figure_3_maxshift_experiment';
imagesc(deltat, angle, mean_eels_exp, [-1 3]);
xlim([-0.5,1.1]);
colormap(my_redblue);
hold on;
for dt = dt_arr
plot(dt * ones(length(0:2:180),1), angle,'LineWidth',2,'LineStyle','-.','Color',[1 1 1]);
end
set_axis_properties(gca,FontSize,FontName,1,0:30:180,-1:0.5:1.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% plot max eels theory
figure;
image_name = 'figure_3_maxshift_theory';
imagesc(t0_vec, angle, mean_eels_sim, [-1 3]);
xlim([-0.5,1.1]);
colormap(my_redblue);
hold on;
for dt = dt_arr
plot(dt * ones(length(0:2:180),1), angle,'LineWidth',2,'LineStyle','-.','Color',[1 1 1]);
end
set_axis_properties(gca,FontSize,FontName,1,0:30:180,-1:0.5:1.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% plot max eels theory experiment comparison
figure;
image_name = 'figure_3_maxshift_comparison';
hold on
ii = 1;
for dt = dt_arr
    [~,ind_0] = min(abs(deltat  - dt));
    error_exp = 0 * angle + .1;
    xx = mean_eels_exp(:,ind_0);
    yy = 0:2:180;
    errorbar(xx(1:3:end),yy(1:3:end),errs_eels_exp(1:3:end, ind_0),'horizontal', ...
      'LineStyle','none','LineWidth',1,Color='black');
    ii = ii + 1;    
end
plot(mean_eels_exp(:,ind_0) * 0,angle,'LineWidth',1,'LineStyle','-.','Color',[1 1 1] * .6);




hold on
ii = 1;
for dt = dt_arr
    [~,ind_0] = min(abs(t0_vec  - dt));
    plot(mean_eels_sim(:,ind_0),0:2:180,LineWidth=3,Color=[0.8 0 0]);
end
set_axis_properties(gca,FontSize,FontName,2,0:30:180,-4:2:4,'','',FontSize,[0 0 0]);
xlim([-5,5]);
ylim([0,180]);
set(gca,"YDir","reverse")
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
box on;
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% spectral density theory

my_yellow = [255, 194, 10] / 255;
my_blue = [12, 123, 220] / 255;
my_orange = [230, 97, 0] / 255;
ind_ang_1 = 5;
ind_ang_2 = 15;
ind_ang_3 = 45;

figure;
image_name = 'spectral_density_theory';
imagesc(f_theory * 1e-12,angle,spectra_theory);
colormap("jet");
xlim([0,2.5]);
hold on;
plot(f_theory * 1e-12, f_theory * 0 + angle(ind_ang_1),'color',my_yellow,'LineStyle','-.','LineWidth',2);
plot(f_theory * 1e-12, f_theory * 0 + angle(ind_ang_2),'color',my_blue,'LineStyle','-.','LineWidth',2);
plot(f_theory * 1e-12, f_theory * 0 + angle(ind_ang_3),'color',my_orange,'LineStyle','-.','LineWidth',2);
set_axis_properties(gca,FontSize,FontName,2,0:30:180,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]);


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% spectral density experiment
figure;
image_name = 'spectral_density_exp';
imagesc(f_exp * 1e-12,angle,spectra_exp);
colormap("jet");
xlim([0,2.5]);
hold on;
plot(f_exp * 1e-12, f_exp * 0 + angle(ind_ang_1),'color',my_yellow,'LineStyle','-.','LineWidth',2);
plot(f_exp * 1e-12, f_exp * 0 + angle(ind_ang_2),'color',my_blue,'LineStyle','-.','LineWidth',2);
plot(f_exp * 1e-12, f_exp * 0 + angle(ind_ang_3),'color',my_orange,'LineStyle','-.','LineWidth',2);
set_axis_properties(gca,FontSize,FontName,2,0:30:180,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

%spectral density comparison
image_name = 'spectral_density_comaprison';


figure; 
hold on
xx = f_exp * 1e-12; 
yy1 =  spectra_exp(ind_ang_1,:);
yy2 =  spectra_exp(ind_ang_2,:);
yy3 =  spectra_exp(ind_ang_3,:);
errors_exp_fourier = 0 * xx + 0.01;
errorbar(xx, yy1, errors_exp_fourier,'vertical','LineStyle','none', ...
    'LineWidth',2,'Color',my_yellow);
errorbar(xx, yy2, errors_exp_fourier,'vertical','LineStyle','none', ...
    'LineWidth',2,'Color',my_blue);
errorbar(xx, yy3, errors_exp_fourier,'vertical','LineStyle','none', ...
    'LineWidth',2,'Color',my_orange);
plot(f_theory * 1e-12, spectra_theory(ind_ang_1,:),'LineWidth',2,'Color',my_yellow)
hold on;
plot(f_theory * 1e-12, spectra_theory(ind_ang_2,:),'LineWidth',2,'Color',my_blue)
plot(f_theory * 1e-12, spectra_theory(ind_ang_3,:),'LineWidth',2,'Color',my_orange)

xlim([0,2.5]);
ylim([0 1])
set_axis_properties(gca,FontSize,FontName,2,0:0.2:1,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
box on;
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

%% color bars for export
% close all
% FontSize = 14;
% image_name = 'colorbar_fourier';
% figure;
% imagesc([0 1])
% colormap("jet")
% colorbar('eastoutside','Ticks',[0 0.5 1]);
% set_axis_properties(gca,FontSize,FontName,2,0:0.2:1,0:0.5:2.5,'','',FontSize,[0 0 0]);
% set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
% exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);
% 
% 
% image_name = 'colorbar_shift';
% figure;
% imagesc([-1 3])
% colormap(my_redblue)
% colorbar('eastoutside','Ticks',[-1 :1: 3]);
% set_axis_properties(gca,FontSize,FontName,2,0:0.2:1,0:0.5:2.5,'','',FontSize,[0 0 0]);
% set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
% exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);
%%
function c = my_redblue(m)
if nargin < 1, m = size(get(gcf,'colormap'),1); end
[r, g, b] = generateRGB(0, 0.25 * m, m, (0:m-1).');
c = [r g b]; 
end

function [r, g, b] = generateRGB(x, y, z, t)
    % Initialize r, g, b as zeros with the same size as t
    r = zeros(size(t));
    g = zeros(size(t));
    b = zeros(size(t));

    % First interval: x <= t <= y
    mask = (t >= x) & (t <= y);
    r(mask) = (t(mask) - x) / (y - x);
    g(mask) = (t(mask) - x) / (y - x);
    b(mask) = 1;

    % Second interval: y < t <= z
    mask = (t > y) & (t <= z);
    r(mask) = 1;
    g(mask) = 1 - ((t(mask) - y) / (z - y));
    b(mask) = 1 - ((t(mask) - y) / (z - y));
end

function [f, spectra] = builder_spectra(psi, energy, time)


    psi_e = psi.*repmat(energy,[length(time),1]);
    psi_exp_e = trapz(energy,psi_e,2)./trapz(energy,psi,2);
    psi_omega = fft(fftshift(psi_exp_e));
    
    pad_psi_exp_e = [psi_exp_e;zeros(5 * length(psi_exp_e),1)];
    dt = (time(2) - time(1))*1e-12;
    fs = 1/dt;
    lsignal = length(pad_psi_exp_e);
    Y = fft(pad_psi_exp_e);
    P2 = abs(Y/lsignal);
    spectra = P2(1:lsignal/2+1);
    f = fs*(0:(lsignal/2))/lsignal;

end

function psi_incoherent = generate_incoherent_spectrum_for_angle(EPD, EOR, TOR, ZOR, eels, w, t_w, e_w, theta)
    
    [t0_vec,eels_comb] = utils_spectrum.calculate_spectrum_from_fields(-cos(2 * theta * pi / 180) * EOR + ...
        EPD, TOR, ZOR);
    psi_incoherent =  eels.incoherent_convolution(utils_spectrum.spectrum_to_coherent_eels(t_w, e_w, eels_comb, t0_vec),...
        w, t_w, e_w);

end


function create_figure_electricfield(T, Z, E, clim, setdir, filename, FontSize)
    figure;
    imagesc(T(:,1), Z(1,:), E, [-clim, clim]);
    set(gca,'FontSize',FontSize);
    xlim([-.3,1.5]);
    ylim([-100,100]);
    xticks(-.3:.3:1.5)
    colormap(utils.redblue);
    pbaspect([2 1 1])
    set(gcf,'position', [200 , 200 , 200 + 150, 200 + 75]);
    colorbar;
    exportgraphics(gcf, [setdir, filename],'resolution', 300);
end



function [TT, ZZ, electric_field_zt,psi_sub_pd ,...
    psi_incoherent_pd,eels, w, e_w, t_w, tt, zz] = electric_field_photodember(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max)

    % set up the eels parameters
    [eels, w, e_w, t_w, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max);
    
    % calculate the photodember interaction voltage
    interact_v_pd_store = utils.eels_pdor(eels,'photodember');
    

    interact_v_pd = circshift(interact_v_pd_store, [0,0]);
    

    alpha_pd_0 =  0.05;
    

    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 ;
    

    [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);

    [TT, ZZ] = ndgrid(tt * 1e12,zz * 1e6);
    

    electric_field_zt = utils_spectrum.derivative_f_dz(interact_v_pd_store, zz, tt);
    electric_field_zt(isnan(electric_field_zt)) = 0;

end

function [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max)

% (TODO: Modify default_parameters)
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

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end


% 
% 
%     [eels_measure, energy,...
%         time] = measurement_plots.data_measurement(angle);
%     
%     deltat = time(time > -1 & time < 1.5);
%     eels_measure = eels_measure(time > -1 & time < 1.5,:);
%     
%     x_c = zeros(1, size(deltat,2));
%     
%     eels_measure = eels_measure ./ sqrt(sum(eels_measure.^2, 2));
%     
%     
%     for i = 1:size(eels_measure,1)
%         [~, ind_max] = max(eels_measure(i,:));
%         x_c(i) = ind_max;
%     end
%     
%     mean_eels = energy(floor(x_c));
% end

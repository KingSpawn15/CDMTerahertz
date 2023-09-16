clear all
close all
setdir = 'figure 3 fourier/results/';


% angle_hwp = angle_hwp_calculator(45);


[deltat, max_eels_measurement] = get_max_eels_measurement(5);
%%
close all
image_name = 'figure_3_maxe';
dt_arr = -0.05;

figure;
FontName = 'ariel';
FontSize = 14;
imagesc(deltat, 0:2:180, max_eels_measurement, [min(max_eels_measurement(:)), max(max_eels_measurement(:))]);
xlim([-0.5,1.5]);
colormap("jet");
hold on;
for dt = dt_arr
plot(dt * ones(length(0:2:180),1), 0:2:180,'LineWidth',2,'LineStyle','-.','Color',[1 1 1]);
end
set_axis_properties(gca,FontSize,FontName,1,0:30:180,-1:0.5:1.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

figure;
image_name = 'figure_3_maxe_section';
hold on
ii = 1;
for dt = dt_arr
    [~,ind_0] = min(abs(deltat  - dt));
    plot(max_eels_measurement(:,ind_0),0:2:180,LineWidth=3,Color='blue');
    leg_str{ii} = num2str(dt);
    ii = ii + 1;
end
plot(max_eels_measurement(:,ind_0) * 0,0:2:180,'LineWidth',1,'LineStyle','-.','Color','black');
set_axis_properties(gca,FontSize,FontName,2,0:30:180,-4:2:4,'','',FontSize,[0 0 0]);
xlim([-5,5]);
ylim([0,180]);
set(gca,"YDir","reverse")
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
box on;
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);




%%
ii = 1;
for ang = 0:2:180

    [max_eels_measure, deltat, energy, eels_measure] = get_measurement_fourier(angle_hwp_calculator(ang));
    
    dataset = eels_measure;
    EnergyCrop = energy;
    time = deltat;
    
    tc = time;
    psi = dataset;
    emat = repmat(EnergyCrop,[length(tc),1]);
    psi_e = psi.*emat;
    psi_exp_e = trapz(EnergyCrop,psi_e,2)./trapz(EnergyCrop,psi,2);
    psi_omega = fft(fftshift(psi_exp_e));
    
    pad_psi_exp_e = [psi_exp_e;zeros(5 * length(psi_exp_e),1)];
    dt = (tc(2) - tc(1))*1e-12;
    fs = 1/dt;
    lsignal = length(pad_psi_exp_e);
    Y = fft(pad_psi_exp_e);
    P2 = abs(Y/lsignal);
    P1 = P2(1:lsignal/2+1);
    P1agg{ii} = P1;
    ii = ii + 1;

end
%%
close all
figure;
image_name = 'spectral_density';
f = fs*(0:(lsignal/2))/lsignal;
fftmat = imagesc(f * 1e-12,0:2:180,filter(cell2mat(P1agg).',4));
colormap("jet");
xlim([0,2.5]);
set_axis_properties(gca,FontSize,FontName,2,0:30:180,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

image_name = 'single_fft';
figure; 
plot(f * 1e-12, P1agg{45},'Marker','+')
xlim([0,2.5]);
ylim([0 1])
set_axis_properties(gca,FontSize,FontName,2,0:0.2:1,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);
%%

angle_90 = 56;

[com_90, deltat, energy, eels_measure] = get_measurement_fourier(angle_90);
    
dataset = eels_measure;
EnergyCrop = energy;
time = deltat;

tc = time;
psi = dataset;
emat = repmat(EnergyCrop,[length(tc),1]);
psi_e = psi.*emat;
psi_exp_e = trapz(EnergyCrop,psi_e,2)./trapz(EnergyCrop,psi,2);
psi_omega = fft(fftshift(psi_exp_e));

dt = (tc(2) - tc(1))*1e-12;
fs = 1/dt;
lsignal = length(tc);
Y = fft(psi_exp_e);
P2 = abs(Y/lsignal);
P1 = P2(1:lsignal/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = fs*(0:(lsignal/2))/lsignal;

%     str_str = strcat('angle',num2str(theta));
data.x = f*1e-12;
data.y = P1;
%     fft_struct.(str_str) = data;
%     ii = ii + 1;

close all;
figure;
set(gcf,'Position',[100,50,600,450]);

plot(f*1e-12,P1,'-.','LineWidth',1,'Color','#6685E1','Marker','+');
xlim([0,5]);
ylim([0,1.2]);
set(gca,'FontSize',14);
xlabel('Frequency (THz)','FontSize',18,'FontName','Times New Roman');
ylabel('Power','FontSize',18,'FontName','Times New Roman');

function [deltat, max_eels_measurement] = get_max_eels_measurement(smoothing)
    
    if nargin < 1
        smoothing = 5;
    end

    ii = 1;
    for angle_pol = 0:2:180
        [max_eels{ii}, deltat, ~, ~] = get_measurement_fourier(angle_hwp_calculator(angle_pol));
        ii = ii + 1;
    end
    
    filter = @(matrix, window)  movmean(movmean(matrix, window, 1) , window, 2);
    max_eels_measurement = filter(cell2mat(max_eels.'), smoothing);

end


function angle_hwp = angle_hwp_calculator(angle_pol, offset)

    if nargin < 2
        offset = 12;
    end
    
    if mod(offset + fix(angle_pol / 2), 2) == 0 
        angle_hwp = offset + fix(angle_pol / 2);
    else
        angle_hwp = offset + fix(angle_pol / 2) - 1;

    end

end



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
    [eels, ~, ~, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj,pd_spot_fwhm, pd_z_max);
    
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

function [com, deltat, energy, eels_measure] = get_measurement_fourier(angle)


    [eels_measure, energy,...
        time] = measurement_plots.data_measurement(angle);
    
    deltat = time;
%     eels_measure = eels_measure;
    
    x_c = zeros(1, size(deltat,2));
    
%     eels_measure = eels_measure ./ sqrt(sum(eels_measure.^2, 2));
    
    
    for i = 1:size(eels_measure,1)
        [~, ind_max] = max(eels_measure(i,:));
        x_c(i) = ind_max;
    end
    
    com = energy(floor(x_c));
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


clear all;
close all;
delete(gcp('nocreate'))
delete(gcp('nocreate'))
parpool(6);




%%
pulse_energy_list = [0.1000    0.1700    0.3000    0.5200    1.0000    1.7300    3.0000    5.1900   10.0000];

i = 1;
for pulse_energy = pulse_energy_list
    [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy);
    interact_v_pd_store_cell{i} = utils.eels_pdor(eels,'photodember');
    i = i+1;
end

interact_v_or_store = utils.eels_pdor(eels,'rectification');


%%

i = 1;
for pulse_energy = pulse_energy_list
    alpha_pd_0 =  .07;
    alpha_or_0 = 20 * 1.3 * 1.2 * pulse_energy / 10;
    
    interact_v_pd = circshift(interact_v_pd_store_cell{i}, [18,0]);
    interact_v_or = circshift(interact_v_or_store, [-15,0]);
    
    t_w = t_w_store-0.2;
    
    
    % alpha_pd = alpha_pd_0; 
    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 + interact_v_or * alpha_or_0;
    [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);
    
    % % alpha_pd = alpha_pd_0; 
    % loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd_0 + interact_v_or * alpha_or_0;
    % [psi_sub_comb , psi_incoherent_comb] = eels.energy_loss_spectrum(loss_spectrum_parameters);
    

    [~, eels_ind_pd] = max(psi_sub_pd,[],2);
    e_exc_pd{i} = e_w(eels_ind_pd);
    i = i+1;
end

%%

[com,deltat,energy,eels_measure] = getmeasurement2();

%%



% [TT, ZZ] = ndgrid(tt,zz);
% vv = interact_v_pd_store;
% electric_field_zt = derivative_f_dz(vv, zz, tt);
% electric_field_zt(isnan(electric_field_zt)) = 0;
% velec = 0.7 * 3 * 10^(8-12);
% [T, Z, ET, t0_vec, eels_cc] = eels_pd(electric_field_zt.', tt * 1e12, zz , velec);
% eels_cc = eels_cc/max(abs(eels_cc(:)))*abs(max(e_exc_pd(:)));
% 
% 
% eels_t = interp1(t0_vec.',eels_cc.',t_w,'linear','extrap');
% psi_assemb_pd  = assemble_psi_sub(t_w, e_w, eels_t, psi_sub_pd);
% psi_incoherent_comb = eels.incoherent_convolution(psi_assemb_pd, w, t_w, e_w);
% 
% 
% 
% 
% % Close all figures
% close all
% 
% % Set directory for saving figures
% setdir = 'Figure 1 fields/results/';

% Set uniform font size and type
FontSize = 14;
FontName = 'helvetica';

% Plot e_exc_pd vs t_w
figure;

plot(e_exc_pd{9}, t_w ,'LineWidth',2)
ylim([-1,1.5]);
xlim([-5,5]);
hold on;
% plot(eels_cc,t0_vec(:),'r--','LineWidth',2)
legend('CDEM','FontSize',FontSize)
hold off
set(gca,'YDir','reverse') ;
set(gca,'FontSize',FontSize);
xticks(-4:2:4);

%%

close all
figure;
hold on;
errors = 0.1 * ones(size(deltat));
colors = turbo(9);
for i = 1:1:9
    plot(e_exc_pd{i} + 2 * (i-1), t_w, 'color', colors(i,:),'LineWidth',1);
    errorbar(com(i,2:3:end) + 2 * (i-1), deltat(2:3:end) - 0.3, 3 * errors(2:3:end), 'horizontal',  'color', colors(i,:),'LineWidth',1);
set(gca, 'YDir', 'reverse');
    hold on;
end
% axis equal;
% daspect([1 1 1]);
ylim([-1, 1.5]);
xlim([-1,20])
ax = gca;
ax.FontSize = 14;
% ax.LineWidth = 2;
ax.YTick = -1:0.5:2;
ax.XTick = [];

stretch_factor = 1.5; % the desired stretching factor
pos = get(gcf, 'Position'); % get the current position vector of the figure
pos(2) = pos(2)*0.5;
pos(3) = pos(3) * stretch_factor; % stretch the height by the stretching factor
pos(4) = pos(4) / stretch_factor;
set(gcf, 'Position', pos); % set the new position vector of the figure
ylabel('\Delta t [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
%%
tiledlayout(1,11,'TileSpacing','compact');
% create a 1 row and 6 column tiled layout for the plots that follow

for i = [1:11]
% loop over the odd-numbered indices from 1 to 11 (i.e., 1, 3, 5, 7, 9, 11)

    if i ~= 1 || i~=11
    % if the current index i is not 1 or 11

        nexttile
        % select the next tile in the layout (i.e., move to the next subplot)

    end
   
%     eels = squeeze(cell2mat(DataSetCropAll(i)))';
    % assign variable eels to the data from DataSetCropAll at index i
    % the data is squeezed and converted to a matrix using cell2mat

    imagesc(energy, deltat, eels_measure{i});
    % create an image plot of the data
    % energy is used for the x-axis, deltat for the y-axis, and eels for the data
    % the transpose of eels is used to orient the plot correctly

    ylim([-1 , 1.5]);
    xlim([-5,5])
    % set limits for the y-axis and x-axis, respectively

    colormap jet
    % set the colormap to 'jet'

    axis square
    % set the aspect ratio of the plot to 1:1

    ax = gca;
%     ax.FontSize = 9;
%     ax.LineWidth = 1;
    ax.YTick = -1:0.5:1.5;
    ax.XTick = -4:2:4;
    % modify the current axis properties, including fontsize, linewidth,
    % y-axis ticks, and x-axis ticks
% 
%     ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
%     xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    % add y-label and x-label to the plot, with specified color and fontsize
    
end



%%
figure;
imagesc(energy, deltat, eels_measure{8});
    % create an image plot of the data
    % energy is used for the x-axis, deltat for the y-axis, and eels for the data
    % the transpose of eels is used to orient the plot correctly

    ylim([-1 , 1.5]);
    xlim([-5,5]);
hold on;
errorbar(com(8,:) + 2 * (1-1), deltat,  2 * errors, 'horizontal',  'color', colors(1,:),'LineWidth',1,'LineStyle','none');

%%

%%%%%%
% % Plot surface of ZZ, TT, and electric_field_zt
% figure;
% s=surf(ZZ, TT, real(electric_field_zt).' / max(real(electric_field_zt(:))));
% s.EdgeColor = 'none';
% 
% % Plot image of tt, zz, and electric_field_zt
% figure;
% imagesc(tt * 1e12, zz * 1e6, real(electric_field_zt)/ max(abs(real(electric_field_zt(:)))),[-1,1]);
% set(gca,'FontSize',FontSize);
% xlim([-.3,1.5]);
% xticks(-.3:.3:1.5)
% colormap(redblue);
% pbaspect([2 1 1])
% colorbar;
% 
% % Plot image of e_w, t_w, and psi_incoherent_comb
% figure;
% imagesc(e_w,t_w, psi_incoherent_comb);
% ylim([-1 , 1.5]);
% colormap jet
% axis square
% set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'\Deltat [ps]','Energy [eV]',FontSize,[0.3 0.3 0.3])
% 
% 
% Plot image of e_w, t_w, and psi_incoherent_pd
%%
figure;
imagesc(e_w,t_w, psi_incoherent_pd);
ylim([-1 , 1.5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'\Deltat [ps]','Energy (eV)',FontSize,[0.3 0.3 0.3])





%% Function to set axis properties
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

function [eels, w, e_w, t_w_store, tt, zz]= eels_setup(pulse_energy_experiment_nj)
%     spot_fwhm = 40e-6;
%     [laser_parameters,discretization_params, utem_parameters,...
%         numerical_parameters] = default_parameters_2(spot_fwhm);

    
    if nargin < 1
        pulse_energy_experiment_nj = 10;
    end

     [laser_parameters,discretization_params, utem_parameters,...
        numerical_parameters] = default_parameters_2();

    laser_parameters.pulse_energy_experiment = 0.1 * pulse_energy_experiment_nj *1e-9;
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

function [com,deltat, energy, eels_measure] = getmeasurement()
    load('saved_matrices\PulseEnergy.mat')
    deltat = Time;
    energy = EnergyCrop;
    % assign variables deltat and energy to values of Time and EnergyCrop, respectively
    % these variables will be used later in the code
    
   
    
    x_c = zeros(11, size(Time,2));
    
    for ii = 1:11
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        eels_measure{ii} = eels;
        % Define the mass of each element in the matrix (assuming unit mass)
        % mass_matrix = ones(size(eels));
        mass_matrix = eels;
        % Initialize a vector to store the x coordinates of the center of mass of each row
       
    
        % Calculate the x coordinate of the center of mass of each row
        for i = 1:size(eels,1)
            row_mass = sum(mass_matrix(i,:));
            for j = 1:size(eels,2)
                x_c(ii,i) = x_c(ii,i) + j*mass_matrix(i,j);
            end
            x_c(ii,i) = x_c(ii,i)/row_mass;
        end
    end

    com = energy(floor(x_c));
end

function [com,deltat, energy, eels_measure] = getmeasurement2()
    load('saved_matrices\PulseEnergy.mat')
    deltat = Time;
    energy = EnergyCrop;
    % assign variables deltat and energy to values of Time and EnergyCrop, respectively
    % these variables will be used later in the code
    
   
    
    x_c = zeros(11, size(Time,2));
    
    for ii = 1:11
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        eels_measure{ii} = eels;
        % Define the mass of each element in the matrix (assuming unit mass)
        % mass_matrix = ones(size(eels));
        mass_matrix = eels;
        % Initialize a vector to store the x coordinates of the center of mass of each row
       
    
        % Calculate the x coordinate of the center of mass of each row
%         for i = 1:size(eels,1)
%             row_mass = sum(mass_matrix(i,:));
%             for j = 1:size(eels,2)
%                 x_c(ii,i) = x_c(ii,i) + j*mass_matrix(i,j);
%             end
%             x_c(ii,i) = x_c(ii,i)/row_mass;
%         end

        for i = 1:size(eels,1)
%             row_mass = sum(mass_matrix(i,:));
            for j = 1:size(eels,2)
                [~, ind_max] = max(mass_matrix(i,:));
                x_c(ii,i) = ind_max;
            end
            x_c(ii,i) = x_c(ii,i);
        end
    end

    com = energy(floor(x_c));
end

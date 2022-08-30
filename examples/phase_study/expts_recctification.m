clearvars;
%%
user_filename  = 'expts_rectification';
dir_saved_matrices = 'examples/phase_study/saved_matrices/';
export_dir = strcat('examples/phase_study/results/',user_filename,'/');
[status, msg, msgID] = mkdir(export_dir);
base_filename = 'eels_';

%%

angle_list = [90];

[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();

laser_parameters.pulse_energy_experiment = 1e-9;
discretization_params.l = 1.5e-12 * 2  * discretization_params.fs;
discretization_params.delay_max = 1.5e-12;

utem_parameters.electron_total_energy = 0.94;
laser_parameters.laser_pulse_time_fwhm = 50e-15;

params.laser_parameters = laser_parameters;
params.discretization_params = discretization_params;
params.utem_parameters =  utem_parameters;
params.numerical_parameters = numerical_parameters;

%%
% build_v_matrices_fn(user_filename, params, angle_list);
%%

%%
% load(strcat(dir_saved_matrices, 'v_struct_',user_filename, '.mat'));
%%
% angle_list = [0, 45, 60, 90, 135, 60, 150];

laser = Laser(laser_parameters);

discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
eels = EELS(eels_parameters);


interaction_gain_factor_rectification = 1;
for interaction_gain_factor_photodember = [0]
    for method =  "rectification"
        for theta_pol_degree = angle_list

            laser.theta_pol =  theta_pol_degree.*(pi/180);
            eels_parameters.laser = laser;
            eels = EELS(eels_parameters);

            loss_spectrum_parameters.method = method;
            loss_spectrum_parameters.interaction_gain_factor_rectification = ...
                interaction_gain_factor_rectification;
            loss_spectrum_parameters.interaction_gain_factor_photodember =...
                interaction_gain_factor_photodember;

            v_struct.(strcat('angle_',num2str(theta_pol_degree))) = ...
                eels.interaction_v(loss_spectrum_parameters);
        end
    end
end

%%
export_dir = 'examples/phase_study/saved_matrices/';
[~,~,~] = mkdir(export_dir);
base_filename = 'v_struct_';
export_filename = strcat(export_dir, base_filename, user_filename,'.mat');
save(export_filename,'v_struct');

%%
for pangle = angle_list


    pol_angle = pangle;
    export_file_name = strcat(export_dir,base_filename,num2str(pol_angle),'.png');

    %%

    discretization = Discretization(discretization_params);
    elec = UTEMElectron(utem_parameters);

    [w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
        discretization.energy, discretization.deltat);


    eels_parameters.numerical_parameters = numerical_parameters;
    laser_parameters.theta_pol = pol_angle * pi / 180;
    laser = Laser(laser_parameters);

    eels_parameters.electron = elec;
    eels_parameters.discretization = discretization;
    eels_parameters.laser = laser;
    eels_parameters.material = IndiumArsenide();
    eels = EELS(eels_parameters);

    interact_v_or = v_struct.(strcat('angle_',num2str(pangle)));
    interact_v_or_store = interact_v_or;
    t_w_store = t_w;

    %%
    alpha_or_0 = 5;

    interact_v_or = circshift(interact_v_or_store, [0,0]);
%     t_w = t_w_store-0.1;


    alpha_pd = 0; alpha_or =  alpha_or_0;
    loss_spectrum_parameters.interact_v = interact_v_or * alpha_or;
    [psi_sub_or , psi_incoherent_or] = eels.energy_loss_spectrum(loss_spectrum_parameters);



    close all;
    tiledlayout('flow');


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


    set(gcf,'Position',[100, 50, 1350, 450*2]);
    exportgraphics(gcf, export_file_name,'resolution' , 400);
end


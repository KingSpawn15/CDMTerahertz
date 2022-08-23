clearvars;
[~ , ~ , ~] = mkdir('examples/phase_study/results/combination/');
load('examples/phase_study/saved_matrices/v_struct_modified.mat');
[status, msg, msgID] = mkdir('examples/phase_study/results');

%%
export_dir = 'examples/phase_study/results/';
base_filename = 'eels_';
pol_angle = 150;
export_file_name = strcat(export_dir,base_filename,num2str(pol_angle),'.png');

%%
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();
laser = Laser(laser_parameters);

discretization = Discretization(discretization_params);
utem_parameters.electron_total_energy = 0.94;
elec = UTEMElectron(utem_parameters);

[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);


eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.numerical_parameters = numerical_parameters;
laser_parameters.laser_pulse_time_fwhm = 500e-15;
laser_parameters.theta_pol = pol_angle * pi / 180;
laser = Laser(laser_parameters);
eels_parameters.laser = laser;
eels_parameters.material = IndiumArsenide();
eels_parameters.material.phase = 0;
eels_parameters.material.gamma_factor = 0.8;
eels = EELS(eels_parameters);

phase = 0;
plot_ind = 1;
close all;
nrow = 2;
ncol = 9;




for delay = [-25, -15, -5, 5, 15, 25]
    close all;
    figure = tiledlayout(nrow,ncol,'Padding', 'none', 'TileSpacing', 'compact');
    kk = 1;
    interaction_gain_factor_photodember =  0.001;
    interaction_gain_factor_rectification = 7*exp(1i*pi/180*0);

    for theta_pol_degree = 10 : 10 : 180

        laser.theta_pol =  theta_pol_degree.*(pi/180);
        eels_parameters.laser = laser;

        eels = EELS(eels_parameters);


        loss_spectrum_parameters.interaction_gain_factor_rectification = ...
            interaction_gain_factor_rectification;
        loss_spectrum_parameters.interaction_gain_factor_photodember =...
            interaction_gain_factor_photodember * exp(1i*phase);
        loss_spectrum_parameters.method = 'combination';
        loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
            v_struct.(strcat('angle_',num2str(theta_pol_degree))) + ...
            interaction_gain_factor_photodember * circshift(v_struct.(strcat('photodember')),[delay 0]);
        [psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);

        nexttile;
        imagesc(e_w,t_w, psi_incoherent);

        ylim([-1 , 1.5])

        if plot_ind == ncol || plot_ind == ncol * nrow
            colorbar;
        end
        colormap jet
        axis square
        ax = gca;
        ax.FontSize = 10;
        ax.LineWidth = 1;
        if (plot_ind == 1 || plot_ind == ncol + 1)
            ax.YTick = -1:0.5:1.5;
            ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',12);
        else
            ax.YTick = [];
        end

        if (plot_ind >= ncol + 1)
            ax.XTick = -4:2:4;
            xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',12);
        else
            ax.XTick = [];
        end


        box on
        plot_ind = plot_ind + 1;

    end

    set(gcf,'Position',[100,100,1600,400]);

    str = ['examples/phase_study/results/combination/','combination_',...
        'pd_gain=',num2str(interaction_gain_factor_photodember),...
        'or_gain=',num2str(interaction_gain_factor_rectification),...
        'delay=',num2str(delay)];
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
end
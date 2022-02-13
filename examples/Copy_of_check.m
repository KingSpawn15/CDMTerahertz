clearvars;
[~ , ~ , ~] = mkdir('results/combination/delay/');
load('saved_matrices/v_struct_5.mat');
[status, msg, msgID] = mkdir('results');

% best parameter 1
% x = [0.0109 0.4259 0]
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();




laser = Laser(laser_parameters);
discretization = Discretization(discretization_params);
elec = UTEMElectron(utem_parameters);

[w, e_w, t_w] = elec.energy_time_grid(numerical_parameters.subsampling_factor,...
    discretization.energy, discretization.deltat);
[fwhm_t , fwhm_e] = utils.calculate_marginals(w, e_w, t_w);

eels_parameters.electron = elec;
eels_parameters.discretization = discretization;
eels_parameters.material = IndiumArsenide_2();
eels_parameters.numerical_parameters = numerical_parameters;



factor_rect = 0;
factor_pd = 7.1e-02; 
delay = 0;

    for interaction_gain_factor_rectification = factor_rect
        for interaction_gain_factor_photodember = factor_pd
            phase = 0;
            plot_ind = 1;
            close all;
            figure = tiledlayout(2,9,'Padding', 'none', 'TileSpacing', 'compact');
            kk = 1;
            for theta_pol_degree = 0
                
                
                laser.theta_pol =  theta_pol_degree.*(pi/180);
                eels_parameters.laser = laser;
                
                eels = EELS(eels_parameters);
                
                
                loss_spectrum_parameters.interaction_gain_factor_rectification = ...
                    interaction_gain_factor_rectification;
                loss_spectrum_parameters.interaction_gain_factor_photodember =...
                    interaction_gain_factor_photodember * exp(1i*phase);
                loss_spectrum_parameters.method = 'photodember';
%                 loss_spectrum_parameters.interact_v = interaction_gain_factor_rectification * ...
%                     v_struct.(strcat('angle_',num2str(theta_pol_degree))) + ...
%                     interaction_gain_factor_photodember * circshift(v_struct.(strcat('photodember')),[delay 0]);
                tic;
                [psi_sub , psi_incoherent] = eels.energy_loss_spectrum(loss_spectrum_parameters);
                save('examples/non essential/phi_incohererent_l.mat','psi_incoherent','psi_sub');
                toc;
                nexttile;
                imagesc(e_w,t_w, psi_incoherent);
                
                ylim([-1 , 1.5]);

                xlim([-3,3]);
                ylim([-0.5,1]);
      
                if plot_ind == 9 || plot_ind == 18
                    colorbar;
                end
                colormap jet
                axis square
                ax = gca;
                ax.FontSize = 10;
                ax.LineWidth = 1;
                if (plot_ind == 1 || plot_ind == 10)
                    ax.YTick = -1:0.5:1.5;
                    ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',12);
                else
                    ax.YTick = [];
                end
                
                if (plot_ind >= 10)
                    ax.XTick = -4:2:4;
                    xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',12);
                else
                    ax.XTick = [];
                end
                
                
                box on
                plot_ind = plot_ind + 1;
                
            end
            set(gcf,'Position',[100,100,1600,400]);
            str = ['results/combination/delay/','combination_',...
                'pd_gain=',num2str(interaction_gain_factor_photodember),...
                'or_gain=',num2str(interaction_gain_factor_rectification)
                ];
%             exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
            
        end
    end






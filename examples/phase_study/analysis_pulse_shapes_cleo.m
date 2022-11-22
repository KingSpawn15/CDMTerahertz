clearvars;
%%
user_filename  = 'addedlongpulse_longsigma=5x_weight=0.2_d';
dir_saved_matrices = 'examples/phase_study/saved_matrices/';
export_dir = strcat('examples/results/cleo_plots/',user_filename,'/');
[status, msg, msgID] = mkdir(export_dir);
base_filename = 'eels_';

%%
[laser_parameters,discretization_params, utem_parameters,...
    numerical_parameters] = default_parameters_2();

laser_parameters.pulse_energy_experiment = 1e-9;
discretization_params.l = 1.5e-12 * 3  * discretization_params.fs;
discretization_params.delay_max = 2 * 1.5e-12;

utem_parameters.electron_total_energy = 0.94;
laser_parameters.laser_pulse_time_fwhm = 650e-15;

params.laser_parameters = laser_parameters;
params.discretization_params = discretization_params;
params.utem_parameters =  utem_parameters;
params.numerical_parameters = numerical_parameters;

%%
% build_v_matrices_fn(user_filename, params);
%%

%%
load(strcat(dir_saved_matrices, 'v_struct_',user_filename, '.mat'));
%%
angle_list = [0, 45, 60, 90, 135, 60, 150];
angle_list = [10 90];

for pangle = angle_list


    pol_angle = pangle;
    export_file_name = strcat(export_dir,base_filename,num2str(pol_angle),'_2.png');

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
    interact_v_pd = v_struct.(strcat('photodember'));
    interact_v_pd_store = interact_v_pd;
    interact_v_or_store = interact_v_or;
    t_w_store = t_w;

    %%

    alpha_pd_0 =  .07;
    alpha_or_0 = 20 * 1.3 * 1.2;

    interact_v_pd = circshift(interact_v_pd_store, [18,0]);
    interact_v_or = circshift(interact_v_or_store, [-15 ,0]);

    %     interact_v_pd = circshift(interact_v_pd_store, [20+8,0]);
    %     interact_v_or = circshift(interact_v_or_store, [-8 +4 ,0]);
    %     t_w = t_w_store-0.2;

    t_w = t_w_store-0.2;
    %     interact_v_or = circshift(v_struct.(strcat('angle_',num2str(angle))),[-15,0]);
    %     interact_v_pd = circshift(v_struct.(strcat('photodember')),[18,0]);


    alpha_pd = alpha_pd_0; alpha_or =  0;
    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
        interact_v_or * alpha_or;
    [psi_sub_pd , psi_incoherent_pd] = eels.energy_loss_spectrum(loss_spectrum_parameters);

    alpha_pd = 0; alpha_or =  alpha_or_0;
    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
        interact_v_or * alpha_or;
    [psi_sub_or , psi_incoherent_or] = eels.energy_loss_spectrum(loss_spectrum_parameters);

    alpha_pd = alpha_pd_0; alpha_or =  alpha_or_0;
    loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
        interact_v_or * alpha_or;
    [psi_sub_com , psi_incoherent_com] = eels.energy_loss_spectrum(loss_spectrum_parameters);

    %%
    close all;
    tiledlayout(1,3,'TileSpacing','compact');



    %     nexttile;
    %     imagesc(e_w,t_w, psi_sub_pd);
    %     ylim([-1 , 1.5]);
    %     colormap jet
    %     axis square
    %     ax = gca;
    %     ax.FontSize = 18;
    %     ax.LineWidth = 1;
    %     ax.YTick = -1:0.2:1.5;
    %     ax.XTick = -4:2:4;
    %     ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    %     xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    %
    %     nexttile;
    %     imagesc(e_w,t_w, psi_sub_or);
    %     ylim([-1 , 1.5]);
    %     colormap jet
    %     axis square
    %     ax = gca;
    %     ax.FontSize = 18;
    %     ax.LineWidth = 1;
    %     ax.YTick = -1:0.2:1.5;
    %     ax.XTick = -4:2:4;
    %     ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    %     xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    %
    %     nexttile;
    %     imagesc(e_w,t_w, psi_sub_com);
    %     ylim([-1 , 1.5]);
    %     colormap jet
    %     axis square
    %     ax = gca;
    %     ax.FontSize = 18;
    %     ax.LineWidth = 1;
    %     ax.YTick = -1:0.2:1.5;
    %     ax.XTick = -4:2:4;
    %     ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    %     xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);

    ax1 = nexttile;
    imagesc(e_w,t_w, psi_incoherent_pd);
    ylim([-1 , 1.5]);
    colormap jet
    axis square
    ax = gca;
    %     ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.YTick = -1:0.5:1.5;
    ax.XTick = -4:2:4;
    %     ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    %     xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    set(ax1.XAxis,'FontName' , 'arial','FontSize',16);
    xlabel('Energy shift (eV)','FontName','arial');
set(ax1.YAxis,'FontName' , 'arial','FontSize',16);
    ax1=  nexttile;
    imagesc(e_w,t_w, psi_incoherent_or);
    ylim([-1 , 1.5]);
    colormap jet
    axis square
    ax = gca;
    %     ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.XTick = -4:2:4;
    ax.YTick = -1:0.5:1.5;
    %     ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    %     xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    set(ax1.XAxis,'FontName' , 'arial','FontSize',16);
    xlabel('Energy shift (eV)','FontName','arial');
    set(ax1.YAxis,'FontName' , 'arial','FontSize',16);
    ax1=  nexttile;
    imagesc(e_w,t_w, psi_incoherent_com);
    ylim([-1 , 1.5]);
    colormap jet
    axis square
    ax = gca;
    %     ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.YTick = -1:0.5:1.5;
    ax.XTick = -4:2:4;
    %     ylabel('\Deltat [ps]','Color',[0.3 0.3 0.3],'FontSize',18);
    %     xlabel('Energy (eV)','Color',[0.3 0.3 0.3],'FontSize',18);
    set(ax1.XAxis,'FontName' , 'arial','FontSize',16);
    set(ax1.YAxis,'FontName' , 'arial','FontSize',16);
    xlabel('Energy shift (eV)','FontName','arial');
    set(gcf,'Position',[100, 50, 900, 300]);
%     exportgraphics(gcf, export_file_name,'resolution' , 400);
end

%%

delay_0 = -15;
delay_90 = -15;
a = (delay_0 + delay_90) /2;
b = (delay_0 - delay_90) /2;

delay_fn = @(theta) fix(a + b * cos(2 * theta));

for delay_or = 100
    nrows = 1;
    ncol = 9;
    close all;
    figure;
    tiledlayout(1,10,'TileSpacing','compact','Padding','compact');
    ii = 1;
    % angle_p_list = sort([[10:10:180],[45,135]]);

    angle_p_list = sort([10 30 50 70 90 110 130 150 170]);


    for angle = angle_p_list


        interact_v_or = v_struct.(strcat('angle_',num2str(angle)));
        interact_v_pd = v_struct.(strcat('photodember'));
        delay_in = delay_fn(angle * pi / 180);

        t_w = t_w_store-0.2;
        interact_v_or = circshift(v_struct.(strcat('angle_',num2str(angle))),[delay_in,0]);
        interact_v_pd = circshift(v_struct.(strcat('photodember')),[18,0]);



        alpha_pd = alpha_pd_0 ; alpha_or =  alpha_or_0*1;
        loss_spectrum_parameters.interact_v = interact_v_pd * alpha_pd + ...
            interact_v_or * alpha_or;
        [psi_sub_com , psi_incoherent_com] = eels.energy_loss_spectrum(loss_spectrum_parameters);

        ax1 = nexttile;
        imagesc(e_w,t_w, psi_incoherent_com);
        pbaspect(ax1,[1 1 1]); %
        colormap('jet');
pbaspect(ax1,[1 1 1]); % <---- move to after-plot
    %     annotation('textbox',[0.15, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
    %         'color',[1 1 1],'LineStyle','none','FontSize',18);
    %     annotation('textbox',[0.15, 0.1, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
    %         'color',[1 1 1],'LineStyle','none','FontSize',8);
    colormap('jet');
    set(ax1.XAxis,'FontName' , 'arial','FontSize',12);
     set(ax1.YAxis,'FontName' , 'arial','FontSize',12);
   
    if ii == 1
        ylabel('\Delta t delay (ps)','FontName','arial');
         
    end

    if ii == 5
       
         xlabel('Energy shift (eV)','FontName','arial');
    end

    %     colorbar;
    %     set(gca,'FontSize',16);
    yticks([]);
    xticks([]);
    if ii > ncol
%         xlabel('Energy shift [eV]','FontName','Times New Roman');
        xticks([-4:2:4]);
    end
xticks([-4:2:4]);
    if (ii == 1 )
%         ylabel('\Delta t (ps)','FontName','Times New Roman');
        yticks([-1 : 0.5 : 1.5]);
    end

    ylim([-1,1.5]);
    xlim([-5,5]);

    %%
    %     str = [dir,...
    %         'single_theta=',num2str(theta),...
    %         ];
    %     exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);

  

        ii = ii + 1;

    end
    set(gcf,'Position',[50,100,1600,200]);
%     set(gcf,'Position',[50,250,1200,450]);
    str = [export_dir,...
        'incoherent_cleo_theory',num2str(delay_or),'_'...
        ];
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);
end
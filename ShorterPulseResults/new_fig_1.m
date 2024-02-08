clear all
close all

fields = optimal_parameters();
EOR = -fields.EOR / 40 * 1.2;
EPD = fields.EPD;
e_w = fields.e_w ;
T = fields.T ;
Z = fields.Z ;
eels_photodember = fields.eels_obj;

%%
close all;
FontName = 'ariel';
FontSize = 15;
clim  = max(abs(EOR(:)));
setdir = 'ShorterPulseResults/results/';
create_figure_electricfield(T, Z, EPD, 2.5e4, setdir, 'field_photodember.png', FontSize);
create_figure_electricfield(T, Z, EOR, clim, setdir, 'field_rectification.png', FontSize);


%%
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

% function [e_w, t_w, psi_incoherent] = generate_incoherent_spectrum_for_angle(EPD, EOR, T, Z, theta, eels_photodember, e_w)
%     
%     [t_w, eels_spectra] = utils_spectrum.calculate_spectrum_from_fields(EPD -cos(2 * theta * pi / 180) * EOR, T, Z*1e-6);
%     psi_assemb = utils_spectrum.spectrum_to_coherent_eels_mod(e_w, t_w, eels_spectra);
%     w = eels_photodember.electron.incoherent_gaussian_blurring_window(e_w,t_w);
%     psi_incoherent =  eels_photodember.incoherent_convolution(psi_assemb, w, t_w, e_w);
% 
% end

function ll = create_line(deltat)
    c = 3*10^(8 - 12 +6);
    ve = 0.7*c;
    z = -100 : 100;
    t = deltat + z / ve;
    ll = line(t,z);
end

function create_figure_electricfield(T, Z, E, clim, setdir, filename, FontSize)
    figure;
    imagesc(T(:,1), Z(1,:), E, [-clim, clim]);
%     h = gcf();
%     h.Renderer = "painters";
    set(gca,'FontSize',FontSize);
    xlim([-.3,1.5]);
    ylim([-100,100]);
    xticks(-.3:.3:1.5)
    colormap(utils.redblue);
    pbaspect([1 1 1])
    set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
    colorbar;
    l0 = create_line(0); l1 = create_line(0.3); l2 = create_line(0.7); 

    set(groot,'defaultAxesXTickLabelRotationMode','manual');
    exportgraphics(gcf, [setdir, filename],'resolution', 300);
%     print(h,'-vector', '-dsvg', [setdir, filename]) 
%     saveas(h, [setdir, filename],'resolution', 300) 
    
end


% function [TPD, ZPD, EPD] = electric_field_photodember(eels_photodember, fitting_parameter_EPD)
%     
%     interact_v_pd = utils.eels_pdor(eels_photodember,'photodember');
%     [TPD, ZPD] = ndgrid(eels_photodember.discretization.t * 1e12,eels_photodember.discretization.z * 1e6);
%     EPD = utils_spectrum.derivative_f_dz(interact_v_pd, eels_photodember.discretization.z, ...
%         eels_photodember.discretization.t) * fitting_parameter_EPD;
%     EPD(isnan(EPD)) = 0;
% 
% end

% function [TOR, ZOR, EOR] = electric_field_rectification(params_rectification, E_max_rectification, delay_or_pd_ps)
%     % close all;
%     % omega = 2 * pi * 1e12;
%     
%     tau = params_rectification.tau ;
%     lambda = params_rectification.lambda ;
%     d = params_rectification.d ;
%     sigma_z = params_rectification.sigma_z ;
%     z = params_rectification.z ;
% 
%     omega_max = 2 * pi * 12e12;
%     np = 10001;
%     
%     materials = opticalresponse;
%     
%     
%     nTHz = @(omega) materials.nTHz_inas_drude(omega / (2 * pi) );
%     refractive_index_data = read_refractive_index('refractive_index_data/InAs.txt');
%     nopt = @(lambda) interpolate_refractive_index(refractive_index_data, lambda * 1e9);
%     delta_lambda = 0.1*1e-9;
%     ngopt = @(lambda) nopt(lambda) - (lambda)*(nopt(lambda + delta_lambda) - ...
%         nopt(lambda))/(delta_lambda);
%     
%     [time_ps, ethz_t, ~, ~] = electric_field_time(lambda, tau, z, d, ngopt, nTHz, nopt, np, omega_max);
%     
%     tt_t = time_ps - delay_or_pd_ps;
% 
% 
%     ethz_t = real(ethz_t(tt_t < 10 & tt_t >-10));
%     tt_t = tt_t(tt_t < 10 & tt_t >-10);
%     vec_z = (-20 : 0.05 : 20 ).'*sigma_z;
% 
%     ethz_t = ethz_t/max(ethz_t);
% %     E_max_rectification = (1.631e6) * 1.65;
%     
%     ethz_t = ethz_t .* E_max_rectification;
%     
%     [TOR, ZOR] = ndgrid(tt_t, vec_z);
%     EOR = repmat(ethz_t, 1, length(vec_z)) .* exp(-ZOR.^2 / (2 * sigma_z^2));
%     EOR = EOR.' ;
%     ZOR = ZOR .* 1e6;
% end
% 
% function [TOR, ZOR, EPDintrap, EOR] = interpolate_field(TOR, ZOR, EOR, TPD, ZPD, EPD)
% 
%     EPDintrap = interp2(TPD.', ZPD.', EPD, TOR.', ZOR.', 'linear', 0);
% 
% end

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
%     axes tight;
end
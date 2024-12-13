clear all
close all

params.e_w = linspace(-5,5,181);
params.exp_theory_time_shift = 0.7;
params.spot_size_fwhm_um_or = 40;

common_fac = .9 * 3e2;
optimal_parameters.weight_pd = 1.2 * 1 * common_fac * 1e20 * (-2);
optimal_parameters.weight_or = 1e4 * .8 * 4 * 2. * common_fac * 50;

% [~, ~, EPD_xz] = get_fields_photodember_meep();
[tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(params.spot_size_fwhm_um_or);
[T, Z] = ndgrid(tc, xc);
e_w = params.e_w;
%%

EOR = ((1/sqrt(3)) * EOR_xz) * optimal_parameters.weight_or /2;
% EOR_ip = - (EOR_xz ) * optimal_parameters.weight_or * (1/sqrt(2));
% EPD = EPD_xz * optimal_parameters.weight_pd;


[t_w_0_200 , psi_incoherent_200] = calculate_incoherent_spectrum_from_fields(EOR, T, Z, e_w, electron_velocity(200));
[t_w_0_150 , psi_incoherent_190] = calculate_incoherent_spectrum_from_fields(EOR, T, Z, e_w, electron_velocity(190));


t_w = t_w_0_200 - params.exp_theory_time_shift + 0.5;
setdir = 'meep_results/results/electron_velocity/';
close all
close all
figure;
image_name = strcat('electron_velocity_comparison');
FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(2,1,"TileSpacing","compact");
ttt.Padding = "loose";

plot_tile(e_w,t_w, psi_incoherent_190);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3]);
plot_tile(e_w,t_w, psi_incoherent_200);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3]);


%%
close all;
FontName = 'ariel';
FontSize = 15;
clim  = max(abs(EOR(:)));
% setdir = 'meep_results/results/electron_velocity';
create_figure_electricfield(T, Z, EOR, clim, setdir, 'field_rectification.png', FontSize);



function ll = create_line(deltat, vel)
    c = 3*10^(8 - 12 +6);
    ve = vel*c;
    z = -120 : 120;
    t = deltat + z / ve;
    ll = line(t,z,'LineWidth',1, 'Color',[0,0.7,0]);
    ll.Color = [ll.Color 0.4];
end

function create_figure_electricfield(T, Z, E, clim, setdir, filename, FontSize, vel)
    figure;
    imagesc(T(:,1), Z(1,:), E, 'Interpolation', 'bilinear',[-clim, clim]);
    set(gca,'FontSize',FontSize);
    xlim([-.3,1.5]);
    ylim([-120,120]);
    yticks(-120:30:120)
    xticks(-.3:.3:1.5)
    colormap(utils.redblue);
    pbaspect([1 1 1])
    set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
    colorbar;
    l0 = create_line(-0.3 + 0.1, electron_velocity(200)); 
    l1 = create_line(0 + 0.1, electron_velocity(200)); 
    l2 = create_line(0.3 + 0.1, electron_velocity(200)); 
    l3 = create_line(-0.3 + 0.1, electron_velocity(100)); 
    set(l3,'Color',[0.7,0,0,0.4])
    l4 = create_line(0 + 0.1, electron_velocity(100)); 
    set(l4,'Color',[0.7,0,0,0.4])
    l5 = create_line(0.3 + 0.1, electron_velocity(100));
    set(l5,'Color',[0.7,0,0,0.4])


    set(groot,'defaultAxesXTickLabelRotationMode','manual');
%     exportgraphics(gcf, [setdir, filename],'resolution', 300);

    
end

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

function v_over_c = electron_velocity(KE_keV)
    % Constants
    m_e = 9.10938356e-31;       % electron mass in kg
    c = 299792458;              % speed of light in m/s
    eV_to_J = 1.60218e-19;      % conversion factor from eV to Joules
    
    % Convert kinetic energy from keV to Joules
    KE = KE_keV * 1e3 * eV_to_J;
    
    % Calculate velocity as a fraction of the speed of light
    v_over_c = sqrt(1 - (m_e * c^2 / (KE + m_e * c^2))^2);
end

function plot_tile(x, y, z)
    FontSize = 18;
    FontName = 'ariel';
    nexttile
    imagesc(x, y, z);
    ylim([-1,1.5]);
    xlim([-5,5]);
    colormap jet
    axis square
    set(groot,'defaultAxesXTickLabelRotationMode','manual');
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3]);
    
end

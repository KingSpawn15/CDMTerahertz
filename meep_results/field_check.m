clear all
close all

fields = optimal_parameters();
EOR = fields.EOR;
EPD = fields.EPD;
e_w = fields.e_w ;
T = fields.T ;
Z = fields.Z ;
eels_photodember = fields.eels_obj;

close all;
FontName = 'ariel';
FontSize = 15;
clim  = max(abs(EOR(:)));
clim_pd  = max(abs(EPD(:)));
setdir = 'meep_results/results/';
create_figure_electricfield(T+.1, Z, EOR, clim, setdir, 'field_rectification_test.png', FontSize);
create_figure_electricfield(T+.1, Z, movmean(movmean(EPD,5,2),2,1), clim_pd, setdir, 'field_photodember_test.png', FontSize);


%%
params.e_w = linspace(-5,5,181);
params.exp_theory_time_shift = 0.7;
params.spot_size_fwhm_um_or = 70;

common_fac = .9 * 3e2;
optimal_parameters.weight_pd = 1.2 * 1 * common_fac ;
optimal_parameters.weight_or = 1e4 * .8 * 4 * 2. * common_fac ;

[~, ~, EPD_xz] = get_fields_photodember_meep();
[tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(params.spot_size_fwhm_um_or);
[T, Z] = ndgrid(tc, xc);
e_w = params.e_w;

[TM, ZM, EPDM] = load_michael_pd();

EOR = EOR_zz * optimal_parameters.weight_or;
EPD = EPD_xz * optimal_parameters.weight_pd;


%%
close all;
FontName = 'ariel';
FontSize = 15;
clim  = max(abs(EPDM(:)));
clim_pd  = max(abs(EPD(:)));
setdir = 'meep_results/results/';
create_figure_electricfield(TM, ZM, EPDM, clim, setdir, 'field_photodember_quasistat.png', FontSize);
create_figure_electricfield(T-.5, Z, movmean(movmean(EPD,5,2),2,1), clim_pd, setdir, 'field_photodember_fdtd.png', FontSize);


%%
% close all
% imagesc(T(:,1)-0.5,Z(1,:),EPD_xz); colormap(utils.redblue); colorbar;
% xlim([-1,1.5]);
% ylim([-100,100]);
% 
% figure;
% imagesc(TM(:,1),ZM(1,:),EPDM); colormap(utils.redblue); colorbar;
% xlim([-1,1.5]);
% ylim([-100,100]);

%%
function [tc, xc, field_pd] = get_fields_photodember_meep()

    load('meep_results\saved_matrices_meep\photodember\combined\field_ez_pd_intensity_10.00t0_0.6fwhm_t_50.mat')
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

end

function [TM, ZM, EPDM] = load_michael_pd()

    load("meep_results\saved_matrices_meep\photodember_michael.mat")
    TM = T;
    ZM = Z;
    EPDM = EPD;
end

function ll = create_line(deltat)
    c = 3*10^(8 - 12 +6);
    ve = 0.7*c;
    z = -100 : 100;
    t = deltat + z / ve;
    ll = line(t,z,'LineWidth',2, 'Color',[0,0.7,0]);
    ll.Color = [ll.Color 0.4];
end

function create_figure_electricfield(T, Z, E, clim, setdir, filename, FontSize)
    figure;
    imagesc(T(:,1), Z(1,:), E, 'Interpolation', 'bilinear',[-clim, clim]);
    set(gca,'FontSize',FontSize);
    xlim([-.3,1.5]);
    ylim([-60,60]);
    yticks(-60:20:60)
    xticks(-.3:.3:1.5)
    colormap(utils.redblue);
    pbaspect([1 1 1])
    set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
    colorbar;
    l0 = create_line(0); l1 = create_line(0.3); l2 = create_line(0.7); 

    set(groot,'defaultAxesXTickLabelRotationMode','manual');
    exportgraphics(gcf, [setdir, filename],'resolution', 300);

    
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

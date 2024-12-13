clear all
close all

% adjust_or = 1.2
% adjust_pd = 1;
% shift = -35;
% spot_size = 120;

params.e_w = linspace(-5, 5, 271);
params.exp_theory_time_shift = 0.1;

common_fac = 1;
optimal_parameters.weight_pd = common_fac * (-0.09 * 1e18);
optimal_parameters.weight_or = common_fac * (1 * 1e9);

[~, ~, EPD_xz] = get_fields_photodember_meep();

e_w = params.e_w;

% Parameter ranges
shift_values = [-32];
adjust_or_values = [1];[1.2];
adjust_pd_values = [0.8];
spotsize_values = [117];


setdir = 'meep_results/results/adjustments/';

for spotsize = spotsize_values
    params.spot_size_fwhm_um_or = spotsize;
    [tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(spotsize);
    [T, Z] = ndgrid(tc, xc);
    for shift = shift_values
        for adjust_or = adjust_or_values
            for adjust_pd = adjust_pd_values
                close all;
                % Apply adjustments
                EPD = shiftWithZeros2D(EPD_xz * optimal_parameters.weight_pd * adjust_pd, 0, shift);
                EOR = (EOR_zz) * optimal_parameters.weight_or * adjust_or;
                EOR_45 = (1/sqrt(2)) * (EOR_xz) * optimal_parameters.weight_or * adjust_or;
                % Calculate spectra
                [t_w_0, psi_incoherent_comb_0] = calculate_incoherent_spectrum_from_fields(-EOR + EPD + EOR_45, T, Z, e_w);
                [~, psi_incoherent_comb_45] = calculate_incoherent_spectrum_from_fields(EPD+ EOR_45, T, Z, e_w);
                [~, psi_incoherent_comb_90] = calculate_incoherent_spectrum_from_fields(EOR + EPD+ EOR_45, T, Z, e_w);
                
                [~, psi_incoherent_or_0] = calculate_incoherent_spectrum_from_fields(-EOR+ EOR_45, T, Z, e_w);
                [~, psi_incoherent_or_45] = calculate_incoherent_spectrum_from_fields(-EOR * 0+ EOR_45, T, Z, e_w);
                [~, psi_incoherent_or_90] = calculate_incoherent_spectrum_from_fields(EOR+ EOR_45, T, Z, e_w);
                
                [~, psi_incoherent_pd] = calculate_incoherent_spectrum_from_fields(EPD, T, Z, e_w);
                
                % Get measurement data
                [~, ~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
                [~, ~, ~, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));
                [~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));
                
                % Time adjustment
                t_w = t_w_0 - params.exp_theory_time_shift;
                
                % Create and save plot
                figure;
                image_name = strcat('tiles_meep_', 'shift', num2str(shift), '_or', num2str(adjust_or), '_pd', num2str(adjust_pd), '_spotsize', num2str(spotsize));
                FontName = 'ariel';
                FontSize = 14;
                ttt = tiledlayout(3, 4, "TileSpacing", "compact");
                ttt.Padding = "loose";
                
                % Create tiles
                plot_tile(energy, deltat, eels_measure_0);
                set_axis_properties(gca, FontSize, FontName, 1, -1:0.5:1.5, [], '', '', FontSize, [0.3 0.3 0.3]);
                plot_tile(e_w, t_w, psi_incoherent_comb_0);
                plot_tile(e_w, t_w, psi_incoherent_pd);
                plot_tile(e_w, t_w, psi_incoherent_or_0);
                
                plot_tile(energy, deltat, eels_measure_45);
                set_axis_properties(gca, FontSize, FontName, 1, -1:0.5:1.5, [], '', '', FontSize, [0.3 0.3 0.3]);
                plot_tile(e_w, t_w, psi_incoherent_comb_45);
                set_axis_properties(gca, FontSize, FontName, 0.01, [], [], '', '', FontSize, [0 0 0, 0]);
                plot_tile(e_w, t_w, psi_incoherent_pd);
                plot_tile(e_w, t_w, psi_incoherent_or_45);
                
                plot_tile(energy, deltat, eels_measure_90);
                set_axis_properties(gca, FontSize, FontName, 1, -1:0.5:1.5, -4:2:4, '', '', FontSize, [0.3 0.3 0.3]);
                plot_tile(e_w, t_w, psi_incoherent_comb_90);
                set_axis_properties(gca, FontSize, FontName, 0.01, [], [], '', '', FontSize, [0 0 0, 0]);
                plot_tile(e_w, t_w, psi_incoherent_pd);
                plot_tile(e_w, t_w, psi_incoherent_or_90);
                
                % Set figure size and export
                set(gcf, 'Position', [200, 200, 200 + 3 * 300, 200 + 400]);
                exportgraphics(gcf, [setdir, image_name, '.png'], 'Resolution', 300);

            end
        end
    end
end

%%

%%
close all;
FontName = 'ariel';
FontSize = 15;
clim  = max(abs(EOR(:)));
clim_pd  = max(abs(EPD(:)));
% setdir = 'meep_results/results/electron_velocity';
create_figure_electricfield(T, Z, EOR, clim, setdir, 'field_rectification.png', FontSize);
create_figure_electricfield(T, Z, EPD, clim_pd, setdir, 'field_photodember.png', FontSize);



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
%     l3 = create_line(-0.3 + 0.1, electron_velocity(100)); 
%     set(l3,'Color',[0.7,0,0,0.4])
%     l4 = create_line(0 + 0.1, electron_velocity(100)); 
%     set(l4,'Color',[0.7,0,0,0.4])
%     l5 = create_line(0.3 + 0.1, electron_velocity(100));
%     set(l5,'Color',[0.7,0,0,0.4])


    set(groot,'defaultAxesXTickLabelRotationMode','manual');
    exportgraphics(gcf, [setdir, filename],'resolution', 300);

    
end

% function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
%     ax.FontSize = FontSize;
%     ax.FontName = FontName;
% %     ax.LineWidth = LineWidth;
%     ax.YTick = YTick;
% 
%     ax.XTick = XTick;
% 
%     ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
%     xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
% %     axes tight;
% end
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

%%
function [tc, xc, field_pd] = get_fields_photodember_meep()

    load('meep_results\saved_matrices_meep\photodember\eps_12\field_ez_pd_intensity_x_10t0_0.5fwhm_t_50sigma_t_40.mat')
%     load('meep_results\saved_matrices_meep\photodember\test\field_ez_pd_intensity_10.00t0_0.5fwhm_t_50.mat')
%     load('meep_results\saved_matrices_meep\photodember\varying_intensities\field_ez_pd_intensity_5.19t0_0.5.mat')
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

end

% function [tc, xc, field_pd] = get_fields_photodember_meep_intensities()
% 
%     load('meep_results/saved_matrices_meep/photodember/spot_size_30_shift04/field_ez_pd_intensity_10t0_0.5.mat')
%     xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
%     tc = tstep : tstep : tstep * size(e_pd,1);
%     field_pd = e_pd.';
% 
% end

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

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end

function shiftedMatrix = shiftWithZeros2D(matrix, shiftRow, shiftCol)
    % Get the size of the matrix
    [rows, cols] = size(matrix);

    % Initialize the output matrix with zeros
    shiftedMatrix = zeros(rows, cols);

    % Calculate the ranges for rows and columns in the shifted matrix
    rowStart = max(1, 1 + shiftRow);
    rowEnd = min(rows, rows + shiftRow);
    colStart = max(1, 1 + shiftCol);
    colEnd = min(cols, cols + shiftCol);

    % Calculate the ranges for rows and columns in the original matrix
    originalRowStart = max(1, 1 - shiftRow);
    originalRowEnd = min(rows, rows - shiftRow);
    originalColStart = max(1, 1 - shiftCol);
    originalColEnd = min(cols, cols - shiftCol);

    % Copy the values from the original matrix to the shifted matrix
    shiftedMatrix(rowStart:rowEnd, colStart:colEnd) = ...
        matrix(originalRowStart:originalRowEnd, originalColStart:originalColEnd);
end

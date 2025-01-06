clear all
close all


params.e_w = linspace(-5, 5, 271);
params.exp_theory_time_shift = 0.1;

common_fac = 1.46 * 0.9281;
optimal_parameters.weight_pd = common_fac * (-0.09 * 1e18);
optimal_parameters.weight_or = common_fac * (1 * 1e9);

[~, ~, EPD_xz] = get_fields_photodember_meep();

e_w = params.e_w;

% Parameter ranges
shift_values = [-32];
adjust_or_values = [1];[1.2];
adjust_pd_values = [0.7];[0.8];
spotsize_values = 117;


setdir = 'meep_results/results/adjustments/';


params.spot_size_fwhm_um_or = spotsize_values;
[tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(spotsize_values);
[T, Z] = ndgrid(tc, xc);


close all;
% Apply adjustments
EPD = shiftWithZeros2D(EPD_xz * optimal_parameters.weight_pd * adjust_pd_values, 0, shift_values);
EOR = (EOR_zz) * optimal_parameters.weight_or * adjust_or_values;



%% Get measurement data

angle = sort([0, 30, 50, 70, 90, 110, 130, 150, 180]);

% [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
for ii = 1:length(angle)
[~, deltat, energy, eels_measure{ii}] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(angle(ii)));
end


%%

for ii = 1:length(angle)
    disp(angle(ii))
    [t_w_0 , psi_incoherent_comb{ii}] = calculate_incoherent_spectrum_from_fields(-EOR * cos(2 * angle(ii) * pi /180) + EPD, T, Z, e_w);
end


%% Get measurement data

[~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
[~, ~, ~, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));
[~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));


%%
    
t_w = t_w_0 - params.exp_theory_time_shift;
setdir = 'meep_results/results/';
figure;
image_name = 'diffusion_model_tiles_multiple_meep';
FontName = 'ariel';
FontSize = 15;
ttt = tiledlayout(2,length(angle),"TileSpacing","compact");
ttt.Padding = "compact";

for ii = 1:length(angle)
plot_tile(energy, deltat, eels_measure{ii});

set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[],'','',FontSize,[0 0 0, 0]);
end

end

for ii = 1:length(angle)
plot_tile(e_w,t_w, psi_incoherent_comb{ii});
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[-4:2:4],'','',FontSize,[0 0 0, 0]);
end
end


set(gcf,'Position',[200,200,200 + 120 * length(angle), 200 +  160]);


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',500);

%%
function [tc, xc, field_pd] = get_fields_photodember_meep()

%     load('meep_results\saved_matrices_meep\photodember\combined\field_ez_pd_intensity_10.00t0_0.6fwhm_t_50.mat')
        load('meep_results\saved_matrices_meep\photodember\eps_12\field_ez_pd_intensity_x_10t0_0.5fwhm_t_50sigma_t_40.mat')

%     load('meep_results\saved_matrices_meep\photodember\test\field_ez_pd_intensity_10.00t0_0.5fwhm_t_50.mat')
%     load('meep_results\saved_matrices_meep\photodember\varying_intensities\field_ez_pd_intensity_5.19t0_0.5.mat')
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

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
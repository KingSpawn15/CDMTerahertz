clear all
close all

%%
[eels_spectra, Energy, Time] = get_measurement_duration();
%%
params.e_w = linspace(-5,5,181);
params.exp_theory_time_shift = 0.6;
params.spot_size_fwhm_um_or = 70;

common_fac = 0.9;
optimal_parameters.weight_pd = 1.2 * 1 * common_fac ;
optimal_parameters.weight_or = 1e4 * .8 * 4 * 2. * common_fac ;

% [~, ~, EPD_xz] = get_fields_photodember_meep();
% [tc, xc, EOR_xz, EOR_yz, EOR_zz] = get_fields_rectification(params.spot_size_fwhm_um_or);


pulse_duration_list = [50, 70, 90, 110, 130, 150, 170, 190, 250, 300, 350];
% pulse_duration_list = [50 150];
for ind = 1 : length(pulse_duration_list)
    filename = strcat('meep_results\saved_matrices_meep\rectification_pulse_time\field_ez',num2str(pulse_duration_list(ind)),'.0_fsshift0.4_ps.mat');
    filename_pd = strcat('meep_results\saved_matrices_meep\photodember\combined\field_ez_pd_intensity_10.00t0_0.6fwhm_t_',num2str(pulse_duration_list(ind)),'.mat');
    [tc, xc, EOR_xz{ind}, EOR_yz, EOR_zz{ind}] = get_fields_rectification_pulseduration(params.spot_size_fwhm_um_or, filename);
    [~, ~, EPD_zz{ind}] = get_fields_photodember_meep_pulse_duration(filename_pd);
end

[T, Z] = ndgrid(tc, xc);
e_w = params.e_w;
%%
for ind = 1 : length(pulse_duration_list)
    EOR{ind} = EOR_zz{ind} * optimal_parameters.weight_or
    EPD{ind} = EPD_zz{ind} * optimal_parameters.weight_pd;
    ETOT{ind} = EOR{ind} + EPD{ind};
    
end
% EOR = EOR_zz * optimal_parameters.weight_or;
% EPD = EPD_xz * optimal_parameters.weight_pd;



%% Get measurement data

% angle = sort([0, 30, 50, 70, 90, 110, 130, 150, 180]);
% 
% % [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
% for ii = 1:length(angle)
% [~, deltat, energy, eels_measure{ii}] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(angle(ii)));
% end


%%

for ii = 1:length(pulse_duration_list)
    disp(pulse_duration_list(ii))
    [t_w_0 , psi_incoherent_comb{ii}] = calculate_incoherent_spectrum_from_fields(ETOT{ii} , T, Z, e_w);
    [t_w_0 , psi_incoherent_or{ii}] = calculate_incoherent_spectrum_from_fields(EOR{ii} , T, Z, e_w);
    [t_w_0 , psi_incoherent_pd{ii}] = calculate_incoherent_spectrum_from_fields(EPD{ii} , T, Z, e_w);
end


%% Get measurement data
% 
% [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
% [~, ~, ~, eels_measure_45] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(45));
% [~, deltat, energy, eels_measure_90] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));


%%
    
t_w = t_w_0 - params.exp_theory_time_shift;
setdir = 'meep_results/results/';
figure;
image_name = 'tiles_multiple_duration';
FontName = 'ariel';
FontSize = 15;
ttt = tiledlayout(4,length(pulse_duration_list),"TileSpacing","compact");
ttt.Padding = "compact";

for ii = 1:length(pulse_duration_list)
plot_tile(Energy, Time - .25, eels_spectra{ii});

set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[],'','',FontSize,[0 0 0, 0]);
end

end

for ii = 1:length(pulse_duration_list)
plot_tile(e_w,t_w, psi_incoherent_comb{ii});
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[-4:2:4],'','',FontSize,[0 0 0, 0]);
end
end

for ii = 1:length(pulse_duration_list)
plot_tile(e_w,t_w, psi_incoherent_or{ii});
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[-4:2:4],'','',FontSize,[0 0 0, 0]);
end
end

for ii = 1:length(pulse_duration_list)
plot_tile(e_w,t_w, psi_incoherent_pd{ii});
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[-4:2:4],'','',FontSize,[0 0 0, 0]);
end
end

set(gcf,'Position',[200,200,200 + 120 * length(pulse_duration_list), 200 +  320]);


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',500);

%%
function [tc, xc, field_pd] = get_fields_photodember_meep_pulse_duration(filename_pd)

    load(filename_pd)
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

end

% function [tc, xc, field_pd] = get_fields_photodember_meep()
% 
%     load('meep_results/saved_matrices_meep/photodember/varying_intensities/field_ez_pd_intensity_10.00t0_0.5.mat')
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

function [eels_spectra, Energy, Time] = get_measurement_duration()


    load("saved_matrices/duration.mat");
    PulseListLength = 11;
    TotalRepetitions = 3;
    Time = -0.75:0.03:3;%[ps]
%     close all
%     figure;
%     FontName = 'ariel';
%     FontSize = 14;
%     ttt = tiledlayout(1,11,"TileSpacing","compact");
%     ttt.Padding = "loose"
    
    for i = 1 : PulseListLength
        
    
        for RepetitionInd = 1:TotalRepetitions          
            
            DataSetAllRep(RepetitionInd,:,:) = squeeze(DataSetCropAll{i,RepetitionInd}(1,:,:));
            Energy = squeeze(EnergyCropAll{i,RepetitionInd}(1,:));
            
        end
        
        dataset = squeeze(mean(DataSetAllRep,1));
        dataset = movmean(dataset,3,2);
        dataset = dataset./trapz(Energy,dataset',2)';
        
        eels_spectra{i} = dataset';
    %     figure(PulseDurationInd)
%         imagesc(Energy,Time,dataset')
    end

end

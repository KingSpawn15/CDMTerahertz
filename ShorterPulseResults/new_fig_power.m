clear all 
close all

delete(gcp('nocreate'));
parpool(6);

% miscellaneous params
e_w = linspace(-5,5,181);

% photodember_parameters
pump_power_nj = 10;
laser_spot_size_fwhm = 40e-6;
fitting_parameter_EPD = (1.26/6.34) * 1.2;

% rectification parameters

params_rectification.tau = 30e-15;
params_rectification.lambda = 800e-9;
params_rectification.d  = .5e-3;
params_rectification.sigma_z = 50e-6;
params_rectification.z = -1e-6;
E_max_rectification = (1.631e6) * 1.65;
delay_or_pd_ps = 0.1;

% pulse energy list 
pulse_energy_list = [0.1000    0.1700    0.3000    0.5200    1.0000    1.7300    3.0000    5.1900   10.0000];

% Rectification field calculation

[TOR, ZOR, EOR_10] = electric_field_rectification(params_rectification, ...
        E_max_rectification, delay_or_pd_ps);

for ii = 1 : length(pulse_energy_list)
    EOR{ii} = EOR_10 * (pulse_energy_list(ii) / 10);
end

% photodember field calculation
for ii = 1 : length(pulse_energy_list)
    eels_photodember = setup_parameters_eels_photodember(pulse_energy_list(ii), laser_spot_size_fwhm);
    [TPD, ZPD, EPD{ii}] = electric_field_photodember(eels_photodember, fitting_parameter_EPD);

end

%%
% photodember field calculation
for ii = 1 : length(pulse_energy_list)
    [T, Z, EPD{ii}, EOR{ii}] = interpolate_field(TOR, ZOR, EOR{ii}, TPD, ZPD, EPD{ii});
end


%% Get measurement data
[com,deltat,energy,eels_measure,errs] = getmeasurement_power();

%%

% [e_w, t_w, psi_incoherent_pd] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, 0 * EOR, T, Z, 45, eels_photodember, e_w);

for ii = 1:length(pulse_energy_list)
    disp(ii)
%     [TPD, ZPD, EPD{ii},psi_sub_pd{ii}, psi_incoherent_pd{ii}, ...
%         eels, w, e_w, t_w, tt, zz] = electric_field_photodember(pulse_energy_list(ii), ...
%         pd_spot_fwhm, pd_z_max);
%     EPD{ii} = EPD{ii}*(1.26/6.34);

[e_w, t_w, psi_incoherent_comb{ii}] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD{ii}, EOR{ii}, T, Z, 90, eels_photodember, e_w);
[t0_vec,eels_comb{ii}] = utils_spectrum.calculate_spectrum_from_fields(EOR{ii} + EPD{ii}, T, Z * 1e-6);

end


%%
close all
figure;
FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(3,3,"TileSpacing","compact");
ttt.Padding = "loose"
nexttile
imagesc(energy, deltat - 0.3, eels_measure{1});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile;
imagesc(energy, deltat- 0.3, eels_measure{8});
    ylim([-1,1.5]);
    xlim([-5,5]);
set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile
imagesc(energy, deltat- 0.3, eels_measure{9});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{1});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3])
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{8});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
axis square
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3])
nexttile
imagesc(e_w,t_w, psi_incoherent_comb{9});
ylim([-1,1.5]);
xlim([-5,5]);
colormap jet
set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3]);
axis square
nexttile([1,3]);
% errors = 0.1 * ones(size(deltat));
colors = turbo(9);
for i = 1:1:9
    plot(eels_comb{i} + 4 * (i-1), t0_vec, 'color', colors(i,:),'LineWidth',1);
    hold on;
    errorbar(com(i,2:2:end) + 4 * (i-1), deltat(2:2:end) - 0.3, ...
       errs(i,2:2:end), 'horizontal',  'color', colors(i,:),'LineStyle','none','LineWidth',.5);
set(gca, 'YDir', 'reverse');
    
end
hold off
set_axis_properties(gca,FontSize,FontName,.2,-1:0.2:1.5,[],'','',FontSize,[0 0 0])
% axis equal;
% daspect([1 1 1]);
ylim([-0.5, 0.5]);
xlim([-4,40]);

set(gcf,'Position',[200,200,200 + 450,200 + 600]); %set paper size (does not affect display)

exportgraphics(gcf, 'new Figures/power_fig.png', 'Resolution',300);

%%
function [com,deltat, energy, eels_measure,errs] = getmeasurement_power()
    load('saved_matrices\PulseEnergy.mat')
    deltat = Time;
    energy = EnergyCrop;

    
    x_c = zeros(11, size(Time,2));
    errs = ones(11, size(Time,2));

    for ii = 1:11
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        eels_measure{ii} = eels;
         
        mass_matrix = eels;
        % Initialize a vector to store the x coordinates of the center of mass of each row
%         errs = zeros(1, size(deltat,2));
        for i = 1:size(eels,1)
%             row_mass = sum(mass_matrix(i,:));
            for j = 1:size(eels,2)
                [~, ind_max] = max(mass_matrix(i,:));
                x_c(ii,i) = ind_max;
            end
            x_c(ii,i) = x_c(ii,i);

            errs(ii,i) = energy(find(mass_matrix(i,:) > 0.7 * max(mass_matrix(i,:)),1,'last')) - ...
            energy(find(mass_matrix(i,:) >  0.7 * max(mass_matrix(i,:)),1,'first'));

        end
    end

    com = energy(floor(x_c));
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
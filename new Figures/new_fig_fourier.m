clear all
close all

fields = optimal_parameters();
EOR = fields.EOR;
EPD = fields.EPD;
e_w = fields.e_w ;
T = fields.T ;
Z = fields.Z ;
eels_photodember = fields.eels_obj;



%% get measurements

smoothing = @(matrix, window)  movmean(movmean(matrix, window, 1) , window, 2);
angle = sort(0:2:180);
mean_eels_measure =  cell(1, length(angle));
eels_measure = cell(1, length(angle));
errs = cell(1, length(angle));
% [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
for ii = 1:length(angle)
[mean_eels_measure{ii}, deltat, energy, eels_measure{ii}, errs{ii}] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(angle(ii)));
end
mean_eels_exp = smoothing(cell2mat(mean_eels_measure.'),4);
errs_eels_exp = cell2mat(errs.');


%% build experiment spectra

spectra_matrix_exp = cell(1, length(angle));
for ii = 1:length(angle)
    [f_exp, spectra_matrix_exp{ii}] = builder_spectra(eels_measure{ii}, energy, deltat(:));
end
spectra_exp = smoothing(cell2mat(spectra_matrix_exp).',4);
spectra_exp = real(spectra_exp) ./ max(abs(real(spectra_exp(:))));

%% base theoretical eels spectra

[~,eels_pd] = utils_spectrum.calculate_spectrum_from_fields(EPD, T, Z * 1e-6);
[t0_vec,eels_or] = utils_spectrum.calculate_spectrum_from_fields(EOR, T, Z * 1e-6);

%% combine eels to get max and spectra
mean_eels_theory = cell(1, length(angle));
psi_incoherent = cell(1, length(angle));
for ii = 1:length(angle)
    disp(angle(ii))
    mean_eels_theory{ii} = -cos(2 * angle(ii) * pi / 180) * eels_or + eels_pd;
    [e_w, t_w, psi_incoherent{ii}] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, EOR, T, Z, angle(ii), eels_photodember, e_w);
end
mean_eels_sim = smoothing(cell2mat(mean_eels_theory).',4);

%%

spectra_matrix_theory= cell(1, length(angle));
for ii = 1:length(angle)
    [f_theory, spectra_matrix_theory{ii}] = builder_spectra(psi_incoherent{ii}, e_w, t_w);
end
spectra_theory = smoothing(cell2mat(spectra_matrix_theory).',4);
spectra_theory = real(spectra_theory) ./ max(abs(real(spectra_theory(:))));


%% plot
close all
setdir = 'new Figures/results/';
FontName = 'ariel';
FontSize = 14;
dt_arr = -0.1;

% plot max eels experiment
figure;
image_name = 'figure_3_maxshift_experiment';
imagesc(deltat, angle, mean_eels_exp, [-1 3]);
xlim([-0.5,1.1]);
colormap(my_redblue);
hold on;
for dt = dt_arr
plot(dt * ones(length(0:2:180),1), angle,'LineWidth',2,'LineStyle','-.','Color',[1 1 1]);
end
set_axis_properties(gca,FontSize,FontName,1,0:30:180,-1:0.5:1.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% plot max eels theory
figure;
image_name = 'figure_3_maxshift_theory';
imagesc(t0_vec, angle, mean_eels_sim, [-1 3]);
xlim([-0.5,1.1]);
colormap(my_redblue);
hold on;
for dt = dt_arr
plot(dt * ones(length(0:2:180),1), angle,'LineWidth',2,'LineStyle','-.','Color',[1 1 1]);
end
set_axis_properties(gca,FontSize,FontName,1,0:30:180,-1:0.5:1.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% plot max eels theory experiment comparison
figure;
image_name = 'figure_3_maxshift_comparison';
hold on
ii = 1;
for dt = dt_arr
    [~,ind_0] = min(abs(deltat  - dt));
    error_exp = 0 * angle + .1;
    xx = mean_eels_exp(:,ind_0);
    yy = 0:2:180;
    errorbar(xx(1:3:end),yy(1:3:end),errs_eels_exp(1:3:end, ind_0),'horizontal', ...
      'LineStyle','none','LineWidth',1,Color='black');
    ii = ii + 1;    
end
plot(mean_eels_exp(:,ind_0) * 0,angle,'LineWidth',1,'LineStyle','-.','Color',[1 1 1] * .6);




hold on
ii = 1;
for dt = dt_arr
    [~,ind_0] = min(abs(t0_vec  - dt));
    plot(mean_eels_sim(:,ind_0),0:2:180,LineWidth=3,Color=[0.8 0 0]);
end
set_axis_properties(gca,FontSize,FontName,2,0:30:180,-4:2:4,'','',FontSize,[0 0 0]);
xlim([-5,5]);
ylim([0,180]);
set(gca,"YDir","reverse")
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
box on;
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% spectral density theory

my_yellow = [255, 194, 10] / 255;
my_blue = [12, 123, 220] / 255;
my_orange = [230, 97, 0] / 255;
ind_ang_1 = 5;
ind_ang_2 = 15;
ind_ang_3 = 45;

figure;
image_name = 'spectral_density_theory';
imagesc(f_theory * 1e-12,angle,spectra_theory);
colormap("jet");
xlim([0,2.5]);
hold on;
plot(f_theory * 1e-12, f_theory * 0 + angle(ind_ang_1),'color',my_yellow,'LineStyle','-.','LineWidth',2);
plot(f_theory * 1e-12, f_theory * 0 + angle(ind_ang_2),'color',my_blue,'LineStyle','-.','LineWidth',2);
plot(f_theory * 1e-12, f_theory * 0 + angle(ind_ang_3),'color',my_orange,'LineStyle','-.','LineWidth',2);
set_axis_properties(gca,FontSize,FontName,2,0:30:180,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]);


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

% spectral density experiment
figure;
image_name = 'spectral_density_exp';
imagesc(f_exp * 1e-12,angle,spectra_exp);
colormap("jet");
xlim([0,2.5]);
hold on;
plot(f_exp * 1e-12, f_exp * 0 + angle(ind_ang_1),'color',my_yellow,'LineStyle','-.','LineWidth',2);
plot(f_exp * 1e-12, f_exp * 0 + angle(ind_ang_2),'color',my_blue,'LineStyle','-.','LineWidth',2);
plot(f_exp * 1e-12, f_exp * 0 + angle(ind_ang_3),'color',my_orange,'LineStyle','-.','LineWidth',2);
set_axis_properties(gca,FontSize,FontName,2,0:30:180,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

%spectral density comparison
image_name = 'spectral_density_comaprison';


figure; 
hold on
xx = f_exp * 1e-12; 
yy1 =  spectra_exp(ind_ang_1,:);
yy2 =  spectra_exp(ind_ang_2,:);
yy3 =  spectra_exp(ind_ang_3,:);
errors_exp_fourier = 0 * xx + 0.01;
errorbar(xx, yy1, errors_exp_fourier,'vertical','LineStyle','none', ...
    'LineWidth',2,'Color',my_yellow);
errorbar(xx, yy2, errors_exp_fourier,'vertical','LineStyle','none', ...
    'LineWidth',2,'Color',my_blue);
errorbar(xx, yy3, errors_exp_fourier,'vertical','LineStyle','none', ...
    'LineWidth',2,'Color',my_orange);
plot(f_theory * 1e-12, spectra_theory(ind_ang_1,:),'LineWidth',2,'Color',my_yellow)
hold on;
plot(f_theory * 1e-12, spectra_theory(ind_ang_2,:),'LineWidth',2,'Color',my_blue)
plot(f_theory * 1e-12, spectra_theory(ind_ang_3,:),'LineWidth',2,'Color',my_orange)

xlim([0,2.5]);
ylim([0 1])
set_axis_properties(gca,FontSize,FontName,2,0:0.2:1,0:0.5:2.5,'','',FontSize,[0 0 0]);
set(gcf,'Position',[200,200,200 + 150,200 +  150]); 
box on;
exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);

%%
function [f, spectra] = builder_spectra(psi, energy, time)


    psi_e = psi.*repmat(energy,[length(time),1]);
    psi_exp_e = trapz(energy,psi_e,2)./trapz(energy,psi,2);
    psi_omega = fft(fftshift(psi_exp_e));
    
    pad_psi_exp_e = [psi_exp_e;zeros(5 * length(psi_exp_e),1)];
    dt = (time(2) - time(1))*1e-12;
    fs = 1/dt;
    lsignal = length(pad_psi_exp_e);
    Y = fft(pad_psi_exp_e);
    P2 = abs(Y/lsignal);
    spectra = P2(1:lsignal/2+1);
    f = fs*(0:(lsignal/2))/lsignal;

end

function c = my_redblue(m)
if nargin < 1, m = size(get(gcf,'colormap'),1); end
[r, g, b] = generateRGB(0, 0.25 * m, m, (0:m-1).');
c = [r g b]; 
end

function [r, g, b] = generateRGB(x, y, z, t)
    % Initialize r, g, b as zeros with the same size as t
    r = zeros(size(t));
    g = zeros(size(t));
    b = zeros(size(t));

    % First interval: x <= t <= y
    mask = (t >= x) & (t <= y);
    r(mask) = (t(mask) - x) / (y - x);
    g(mask) = (t(mask) - x) / (y - x);
    b(mask) = 1;

    % Second interval: y < t <= z
    mask = (t > y) & (t <= z);
    r(mask) = 1;
    g(mask) = 1 - ((t(mask) - y) / (z - y));
    b(mask) = 1 - ((t(mask) - y) / (z - y));
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
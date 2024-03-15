clear all
close all

e_w = linspace(-5,5,181);
[tc, xc, EPD] = photodember_meep();

% [t_w, eels_spectra] = utils_spectrum.calculate_spectrum_from_fields(EPD -cos(2 * theta * pi / 180) * EOR, T, Z*1e-6);

[TPD, ZPD] = ndgrid(tc, xc);

%%

EPD_smooth = movmean(movmean(EPD,3,1),3,2);





[t_w, eels_spectra] = utils_spectrum.calculate_spectrum_from_fields(EPD, TPD, ZPD*1e-6);
[~, eels_spectra_smooth] = utils_spectrum.calculate_spectrum_from_fields(EPD_smooth, TPD, ZPD*1e-6);

psi_assemb = utils_spectrum.spectrum_to_coherent_eels_mod(e_w, t_w, eels_spectra_smooth * 1e2 * 3);


%%

eels_photodember = setup_parameters_eels_photodember(10, 40e-6);
w = eels_photodember.electron.incoherent_gaussian_blurring_window(e_w,t_w);
psi_incoherent =  eels_photodember.incoherent_convolution(psi_assemb, w, t_w, e_w);

%%
close all
figure;

FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(2,3,"TileSpacing","compact");
ttt.Padding = "loose";

nexttile
imagesc(TPD(:,1),ZPD(1,:),EPD);
colormap(utils.redblue);
axis square
colorbar;

nexttile;
imagesc(TPD(:,1),ZPD(1,:),EPD_smooth);
colormap(utils.redblue);
axis square
colorbar;

nexttile;
plot(eels_spectra, t_w - .3);
set(gca,'YDir','reverse');
ylim([-1,1.5]);
xlim([-4 * max(abs(eels_spectra)), 4 * max(abs(eels_spectra))])
axis square

nexttile;
plot(eels_spectra_smooth, t_w -.3);
axis square;
ylim([-1,1.5]);
xlim([-4 * max(abs(eels_spectra)), 4 * max(abs(eels_spectra))])
set(gca,'YDir','reverse');

nexttile;
imagesc(e_w, t_w - .3, psi_assemb);
colormap("jet")
axis square;
ylim([-1,1.5]);

nexttile;
imagesc(e_w, t_w - .3, psi_incoherent);
colormap("jet")
axis square;
ylim([-.5,1]);
xlim([-3,3]);
% xlim([-4 * max(abs(eels_spectra)), 4 * max(abs(eels_spectra))])
% set(gca,'YDir','reverse');

set(gcf,'Position',[200,200,200 + 3 * 300,200 +  400]);


%%
    


function [tc, xc, field_pd] = photodember_meep()
    load('article_check/saved_matrices_meep/photodember/spot_size_30/field_ez_pd_intensity_5.mat')
    

  
    xc =  - zstep * size(e_pd,2) / 2 : zstep:  zstep * size(e_pd,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_pd,1);
    field_pd = e_pd.';

end

results_dir = 'InAs_properties_studies/results/';
materials = opticalresponse;
nTHz = @(omega) materials.nTHz_inas_drude(omega / (2 * pi) );

omega_min_thz = 0;
omega_max_thz = 200;

omega = linspace(omega_min_thz, omega_max_thz, 100) * 1e12;

n_real = real(nTHz(omega));
n_imag = imag(nTHz(omega));

%%
close all;
figure;
tiledlayout(1,2)
nexttile;
plot(omega/1e12/2/pi, n_real.^2 - n_imag.^2 );
hold on;
plot(omega/1e12/2/pi, omega * 0, 'k--' );
ylim([-200,100]);
nexttile;
plot(omega/1e12/2/pi, n_imag);

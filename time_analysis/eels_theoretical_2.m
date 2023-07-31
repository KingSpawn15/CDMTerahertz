function [t0_vec, eels] = eels_theoretical_2(tau, lambda, d, z, sigma_z)
% close all;
% omega = 2 * pi * 1e12;
omega_max = 2 * pi * 12e12;
np = 10001;

materials = opticalresponse;


nTHz = @(omega) materials.nTHz_inas_drude(omega / (2 * pi) );
refractive_index_data = read_refractive_index('refractive_index_data/InAs.txt');
nopt = @(lambda) interpolate_refractive_index(refractive_index_data, lambda * 1e9);
delta_lambda = 0.1*1e-9;
ngopt = @(lambda) nopt(lambda) - (lambda)*(nopt(lambda + delta_lambda) - ...
    nopt(lambda))/(delta_lambda);

[time_ps, ethz_t, ~, ~] = electric_field_time(lambda, tau, z, d, ngopt, nTHz, nopt, np, omega_max);

velec = 0.7 * 3 * 10^(8-12);

tt_t = time_ps;


ethz_t = real(ethz_t(tt_t < 10 & tt_t >-10));
ethz_t = ethz_t/max(ethz_t);
tt_t = tt_t(tt_t < 10 & tt_t >-10);
[~, ~, ~, t0_vec, eels] = eels_calc(ethz_t, tt_t, sigma_z, velec);
eels = -eels./max(abs(eels));


end



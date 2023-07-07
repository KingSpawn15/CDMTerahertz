clear all;
refractive_index_data = read_refractive_index('refractive_index_data/InAs.txt');
% [n,k] = interpolate_refractive_index(refractive_index_data, 500)

delta_lambda = 0.1;
nopt = @(lambda) interpolate_refractive_index(refractive_index_data, lambda);
alpha_opt = @(lambda) interpolate_absorption(refractive_index_data, lambda);
ngopt = @(lambda) real(nopt(lambda)) - (lambda)*(real(nopt(lambda + delta_lambda)) - ...
    real(nopt(lambda)))/delta_lambda;

%%



no_800 = nopt(800);
ng_800 = ngopt(800);
alpha_800 = alpha_opt(800) * 1e9 * 1e-2;


materials = opticalresponse;
nTHz = @(omega) materials.nTHz_inas(omega / (2 * pi) );

nu = (0:0.1:10).' * 1e12;
plot(nu, real(sqrt(nTHz(2 * pi * nu))), nu, imag(sqrt(nTHz(2 * pi * nu))))
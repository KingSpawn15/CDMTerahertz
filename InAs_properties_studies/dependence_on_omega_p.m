clear all
close all

lambda = 800e-9;
tau = 30e-15;
z = -1e-6;
d  = .5e-3;

materials = opticalresponse;
refractive_index_data = read_refractive_index('refractive_index_data/InAs.txt');
nopt = @(lambda) interpolate_refractive_index(refractive_index_data, lambda * 1e9);
delta_lambda = 0.1*1e-9;
ngopt = @(lambda) nopt(lambda) - (lambda)*(nopt(lambda + delta_lambda) - ...
    nopt(lambda))/(delta_lambda);


np = 10001;
omega_max = 2 * pi * 12e12;

i = 1;

time_ps_c = {};
ethz_t_c = {};
for factor = 0.2:0.5:10
    disp(i)
    nTHz = @(omega) materials.nTHz_inas_drude_modified(omega / (2 * pi) , (78.1/14.7));
    [time_ps, ethz_t, omg, eomg] = electric_field_time(lambda, tau, z, d, ...
        ngopt, nTHz, nopt, ...
        np, omega_max);
    ethz_t_re = real(ethz_t);
    time_ps_c{i} = time_ps;
    ethz_t_c{i} = ethz_t_re;
    i = i + 1;
end

%%
close all
plot(cell2mat(time_ps_c), cell2mat(ethz_t_c))
xlim([-2,10]);
set(gca,'TickLength',[.02,.02]);



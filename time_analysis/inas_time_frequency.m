clear all
close all
tau = 150e-15;
omega =  (0 : 2 * pi * 0.05e12 : 2 * pi * 20e12).';
z = -1e-6;
d_arr = [.05, .5, 1, 2, 5] * 1e-3;

close all;
omega_max = 2 * pi * (12)*1e12;
np = 10001;

materials = opticalresponse;


nTHz = @(omega) materials.nTHz_inas_drude(omega / (2 * pi) );


refractive_index_data = read_refractive_index('refractive_index_data/InAs.txt');
nopt = @(lambda) (interpolate_refractive_index(refractive_index_data, lambda * 1e9));
delta_lambda = 0.1*1e-9;
ngopt = @(lambda) nopt(lambda) - (lambda)*(nopt(lambda + delta_lambda) - ...
    nopt(lambda))/(delta_lambda);


%%
lambda_arr = [740, 770, 800, 860, 1560]*1e-9;

alpha1 = @(lambda)  (2 * pi / lambda) * imag(nopt(lambda)) * 1e1;
alpha2 = @(lambda)  (2 * pi / lambda) * imag(nopt(lambda)) * 1e0;
alpha3 = @(lambda)  (2 * pi / lambda) * imag(nopt(lambda)) * 1e-1;
alpha4 = @(lambda)  (2 * pi / lambda) * imag(nopt(lambda)) * 1e-2;
alpha5 = @(lambda)  (2 * pi / lambda) * imag(nopt(lambda)) * 1e-3;

tau_arr = [200, 150, 100, 50, 30]*1e-15;
[time_ps1, ethz_t1, omg1, eomg1] = electric_field_time(lambda_arr(1), tau_arr(1), ...
    z, d_arr(2), ngopt, nTHz, nopt, np, omega_max);
[time_ps2, ethz_t2, omg2, eomg2] = electric_field_time(lambda_arr(2), tau_arr(1), ...
    z, d_arr(2), ngopt, nTHz, nopt, np, omega_max);
[time_ps3, ethz_t3, omg3, eomg3] = electric_field_time(lambda_arr(3), tau_arr(1),...
    z, d_arr(2), ngopt, nTHz, nopt, np, omega_max);
[time_ps4, ethz_t4, omg4, eomg4] = electric_field_time(lambda_arr(4), tau_arr(1),...
    z, d_arr(2), ngopt, nTHz, nopt, np, omega_max);
[time_ps5, ethz_t5, omg5, eomg5] = electric_field_time(lambda_arr(5), tau_arr(1),...
    z, d_arr(2), ngopt, nTHz, nopt, np, omega_max);

plot_scaling = 0.3;

figure;

plot(time_ps1, real(ethz_t1)/max(real(ethz_t1))*plot_scaling + 0.5 * 4 , ...
    time_ps2, real(ethz_t2)/max(real(ethz_t2))*plot_scaling+ 0.5 * 3,...
     time_ps3, real(ethz_t3)/max(real(ethz_t3))*plot_scaling+ 0.5 * 2,...
      time_ps4, real(ethz_t4)/max(real(ethz_t4))*plot_scaling+ 0.5 * 1,...
       time_ps5, real(ethz_t5)/max(real(ethz_t5))*plot_scaling+ 0.5 * 0,...
    'LineWidth',1);
hold on;
expo = @(t,tau) exp(-(t.^2)/(tau^2)/2);
tt_expo = (-2:0.01:3).';
plot(tt_expo, expo(tt_expo, tau_arr(1) * 1e12)*plot_scaling + 0.5 * 4 , ...
    tt_expo, expo(tt_expo, tau_arr(1) * 1e12)*plot_scaling + 0.5 * 3 , ...
    tt_expo, expo(tt_expo, tau_arr(1) * 1e12)*plot_scaling + 0.5 * 2 , ...
    tt_expo, expo(tt_expo, tau_arr(1) * 1e12)*plot_scaling + 0.5 * 1 , ...
    tt_expo, expo(tt_expo, tau_arr(1) * 1e12)*plot_scaling + 0.5 * 0 , ...
    'LineWidth',1,'LineStyle','--','Color',[0.3,0.3,0.3]);
hold off
xlim([-1,2.5]);


exportgraphics(gcf,'results_frequency/inas_time_frequency.png','Resolution',500)

%%
figure;

plot_scaling = 0.3;

plot(omg1/(2 * pi * 1e12), abs(eomg1).^2/max(abs(eomg1).^2)*plot_scaling + 0.5 * 4 , ...
    omg2/(2 * pi * 1e12), abs(eomg2).^2/max(abs(eomg2).^2)*plot_scaling + 0.5 * 3 , ...
    omg3/(2 * pi * 1e12), abs(eomg3).^2/max(abs(eomg3).^2)*plot_scaling + 0.5 * 2 , ...
    omg4/(2 * pi * 1e12), abs(eomg4).^2/max(abs(eomg4).^2)*plot_scaling + 0.5 * 1 , ...
    omg5/(2 * pi * 1e12), abs(eomg5).^2/max(abs(eomg5).^2)*plot_scaling + 0.5 * 0 , ...
    'LineWidth',1);

hold on;

expo_omg = @(omg,tau) exp(-(omg.^2 .* tau^2)/2);
nu_expo = (0:0.01:10).';
plot(nu_expo,expo_omg(nu_expo *(2 * pi * 1e12), tau_arr(5))*plot_scaling + 0.5 * 4 , ...
    nu_expo, expo_omg(nu_expo *(2 * pi * 1e12), tau_arr(5))*plot_scaling + 0.5 * 3 , ...
    nu_expo , expo_omg(nu_expo *(2 * pi * 1e12), tau_arr(5))*plot_scaling + 0.5 * 2 , ...
    nu_expo, expo_omg(nu_expo *(2 * pi * 1e12), tau_arr(5))*plot_scaling + 0.5 * 1 , ...
    nu_expo, expo_omg(nu_expo *(2 * pi * 1e12), tau_arr(5))*plot_scaling + 0.5 * 0 , ...
    'LineWidth',1,'LineStyle','--','Color',[0.3,0.3,0.3]);
hold off
xlim([0,12])
ylim([-0.5,2.5])
pbaspect([1 3 1])
exportgraphics(gcf,'results_frequency/inas_omega_frequency.png','Resolution',500)

%%

%%
% close all
velec = 0.7 * 3 * 10^(8-12);
sigma_z = 55 * 1e-6;

ethz_t = ethz_t5;
tt_t = time_ps5;

% ethz_t = real(ethz_t(tt_t < 3 & tt_t >-2));
% tt_t = tt_t(tt_t < 3 & tt_t >-2);

 
ethz_t = real(ethz_t(tt_t < 10 & tt_t >-10));
ethz_t = ethz_t/max(ethz_t);
tt_t = tt_t(tt_t < 10 & tt_t >-10);
[T, Z, ET, t0_vec, eels] = eels_calc(ethz_t, tt_t, sigma_z, velec);

ET = 5 * ET;
% 
figure;

s = surf(T, Z, ET);
s.EdgeColor = 'none';

% s.EdgeAlpha = 0.99;
% s.FaceAlpha = 0.9;
brighten(-.5);
hold on;

eels = real(eels);
c = line((t0_vec), (t0_vec) * 0 + 5*sigma_z, 4+2.5 * eels / max(abs(eels)));
set(c,'Color','blue')
c.LineWidth = 1;
t_grid = 0;
[line_vert3, line_horz3] = grid_lines(t_grid, 5 * sigma_z, t0_vec, eels);
[line_electron_0] = electron_travel(t_grid, 3 * sigma_z, T, Z, ET,velec);
shift = 9;
tau_30 = 30 * 10^(-15 + 12);
tvec = (-0.5:0.002:0.5).'; 
[i_line] = pulse(tau_30, tvec, 5 * sigma_z, shift);



xlim([-1,1.5])
ylim([-5,6]*sigma_z)
zticks([])
view([-4,19])
hold off;

exportgraphics(gcf, "results_frequency/loss_spec_inas_frequency.png",'Resolution',500)

function [line_electron] = electron_travel(t0, zmax, T, Z, ET, velec)
%     t0 = -2;
    zz_l = linspace(-zmax , zmax, 50).';
    tt_l = zz_l / velec + t0;
    et_l = interp2(T', Z', ET', tt_l, zz_l);
    line_electron = line(tt_l, zz_l, et_l);
    set(line_electron,'Color','red');
    line_electron.LineWidth = 2;

end

function [i_line] = pulse(tau, tvec, max_z, shift)

    if nargin < 4
        shift = 2;
    end
    
    it = 2 * exp(-(tvec.^2) ./ tau.^2 /2) + shift;
    i_line = line(tvec, tvec.*0 + max_z, it);
    i_line.Color = 'red';
    
    i_line.LineWidth = 1;

end

function [line_vert, line_horz] = grid_lines(tp, max_z, t0_vec, eels)
    
    eels_plot_fn = @(tp) interp1(t0_vec,4-2.5 * eels / max(abs(eels)), tp);
    
    
    eels_fn_0 = eels_plot_fn(tp);
    zp = (0: 0.1 * eels_fn_0 : eels_fn_0).';
    xp = zp * 0 + tp;
    yp = zp * 0 + max_z;
    yp_h = linspace(0, max_z , length(zp)).';
    zp_h = yp_h * 0;
    
    
    line_vert = line(xp, yp, zp);
    line_vert.Color = 'black';
    line_vert.LineStyle = '--';
    line_vert.LineWidth = 1;
    
    line_horz = line(xp, yp_h, zp_h);
    line_horz.Color = 'black';
    line_horz.LineStyle = '--';
    line_horz.LineWidth = 1;

end




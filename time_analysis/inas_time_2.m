clear all
close all
tau = 150e-15;
omega =  (0 : 2 * pi * 0.05e12 : 2 * pi * 10e12).';
z = -1e-6;
d = .5e-3;

close all;
omega_max = 2 * pi * 3e12;
np = 10001;

materials = opticalresponse;


nTHz = @(omega) (materials.nTHz_inas(omega / (2 * pi) ));
refractive_index_data = read_refractive_index('refractive_index_data/InAs.txt');
nopt = @(lambda) (interpolate_refractive_index(refractive_index_data, lambda * 1e9));
delta_lambda = 0.1*1e-9;
ngopt = @(lambda) nopt(lambda) - (lambda)*(nopt(lambda + delta_lambda) - ...
    nopt(lambda))/(delta_lambda);

% nTHz = @(omega) materials.nTHz_znte(omega / (2 * pi) * 1e-12);
% nopt = @(lambda) materials.nopt_znte(1e6 * lambda);
% ngopt = @(lambda) materials.ngopt_znte(lambda * 1e6);

lambda_arr = [740, 770, 800, 860, 1560]*1e-9;
% [time_ps1, ethz_t1] = electric_field_time(lambda_arr(1), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
% [time_ps2, ethz_t2] = electric_field_time(lambda_arr(2), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
% [time_ps3, ethz_t3] = electric_field_time(lambda_arr(3), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
% [time_ps4, ethz_t4] = electric_field_time(lambda_arr(4), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
% [time_ps5, ethz_t5] = electric_field_time(lambda_arr(5), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
% % 
% plot_scaling = 0.25;
% 
% plot(time_ps1, real(ethz_t1)/max(real(ethz_t1))*plot_scaling + 0.5 * 4 , ...
%     time_ps2, real(ethz_t2)/max(real(ethz_t1))*plot_scaling+ 0.5 * 3,...
%      time_ps3, real(ethz_t3)/max(real(ethz_t1))*plot_scaling+ 0.5 * 2,...
%       time_ps4, real(ethz_t4)/max(real(ethz_t1))*plot_scaling+ 0.5 * 1,...
%        time_ps5, real(ethz_t5)/max(real(ethz_t1))*plot_scaling+ 0.5 * 0,...
%     'LineWidth',2);
% xlim([-2,3])


% figure;
% 
tau_arr = [30, 150, 100, 50, 30]*1e-15;
[time_ps1, ethz_t1] = electric_field_time(lambda_arr(3), tau_arr(1), z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps2, ethz_t2] = electric_field_time(lambda_arr(3), tau_arr(2), z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps3, ethz_t3] = electric_field_time(lambda_arr(3), tau_arr(3), z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps4, ethz_t4] = electric_field_time(lambda_arr(3), tau_arr(4), z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps5, ethz_t5] = electric_field_time(lambda_arr(3), tau_arr(5), z, d, ngopt, nTHz, nopt, np, omega_max);
% 
plot_scaling = 0.3;

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
    tt_expo, expo(tt_expo, tau_arr(2) * 1e12)*plot_scaling + 0.5 * 3 , ...
    tt_expo, expo(tt_expo, tau_arr(3) * 1e12)*plot_scaling + 0.5 * 2 , ...
    tt_expo, expo(tt_expo, tau_arr(4) * 1e12)*plot_scaling + 0.5 * 1 , ...
    tt_expo, expo(tt_expo, tau_arr(5) * 1e12)*plot_scaling + 0.5 * 0 , ...
    'LineWidth',1,'LineStyle','--','Color',[0.3,0.3,0.3]);



hold off
xlim([-2,3])
exportgraphics(gcf,'inas_time.png','Resolution',500)
%%

%%
% close all
velec = 0.7 * 3 * 10^(8-12);
sigma_z = 20 * 1e-6;

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


c = line((t0_vec), (t0_vec) * 0 + 5*sigma_z, 4-2.5 * eels / max(abs(eels)));
set(c,'Color','blue')
c.LineWidth = 1;

% [line_vert1, line_horz1] = grid_lines(-2, 5 * sigma_z, t0_vec, eels);
% [line_vert2, line_horz2] = grid_lines(-1.5, 5 * sigma_z, t0_vec, eels);
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

exportgraphics(gcf, "loss_spec_inas.png",'Resolution',500)
% figure;
% plot(2.5 * -eels / max(abs(eels)), (t0_vec) ,'Marker','+')
% ylim([-1,1.5])
% xlim([-4.5,4.5])
% set(gca, 'YDir','reverse')
%%
% plot((t0_vec), -eels / max(abs(eels)))
% xlim([-1,1.5])
%%
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

function [time_ps, ethz_t] = electric_field_time(lambda, tau, z, d, ngopt, nTHz, nopt, np, omega_max)
    
    if nargin < 8
        omega_max = 2 * pi * 2e12;
        if nargin < 7
            np = 1001;
        end
    end


    C = 3e8;
    mu0 = 4 * pi * 1e-7;
    
    
    Ii = @(omega, tau) tau * sqrt(pi) * exp(-(tau^2 * omega.^2)/2);
    
    q = @(omega) omega .* nTHz(omega) / C;
    q0 = @(omega, lambda) omega .* ngopt(lambda) / C;
    qv = @(omega) omega / C;
    
    
    
    S = @(omega, lambda, tau) (1i * mu0 * omega.^2 .* Ii(omega, tau)) ./ (C * nopt(lambda) .* (q(omega) + ...
        q0(omega, lambda)));
    
    r = @(omega) (q(omega) - qv(omega)) ./ (q(omega) + qv(omega));
    t = @(omega) (2 * q(omega)) ./ (q(omega) + qv(omega));
    L = @(omega, lambda, z, d) (1i * (-exp(1i * z * q(omega)) + exp(1i * d * q0(omega, lambda)))) ./ (q(omega) - ...
        q0(omega, lambda));
    emin = @(omega, lambda, tau, z, d) (exp(1i*(d - z)*q(omega)) .* (L(omega, lambda, d, d) .* r(omega) ...
        + (0.5i * (1 + r(omega)) .* (-exp(1i * d * q0(omega, lambda)) + ...
        exp(1i * d * q(omega)) .* r(omega))) ./ q(omega)) .* S(omega, lambda, tau)) ./ (1 - ...
        exp(2 * 1i * d * q(omega)) .* r(omega).^2);
    
    epl = @(omega, lambda, tau, z, d) L(omega, lambda, z, d) .* S(omega, lambda, tau) + ...
        exp(1i * z * q(omega)) .* r(omega) .* ((0.5i * (1 + ...
        r(omega)) .* S(omega, lambda, tau)) ./ q(omega) + ...
        (exp(1i * d * q(omega)) .* (L(omega, lambda, d, d) .* r(omega) + (0.5i * (1 + ...
        r(omega)) .* (-exp(1i * d * q0(omega, lambda)) + ...
        exp(1i * d * q(omega)) .* r(omega))) ./ q(omega)) .* S(omega, lambda, tau)) ./ (1 - ...
        exp(2 * 1i * d * q(omega)) .* r(omega).^2));
    
    ef = @(omega, lambda, tau, z, d) epl(omega, lambda, tau, d, d) .* t(omega) - (0.5i * exp(1i * d * q0(omega, lambda) ...
        - 1i * (-d + z) * qv(omega)) .* S(omega, lambda, tau) .* t(omega)) ./ q(omega);
    
%     eb = @(omega, lambda, tau, z, d) emin(omega, lambda, tau, 0, d) .* t(omega) + ...
%         (0.5i * S(omega, lambda, tau) .* t(omega)) ./ (exp(1i * z * qv(omega)) .* q(omega));
    
     eb = @(omega, lambda, tau, z, d) exp(1i*z*qv(omega)).*(emin(omega,lambda,tau,0,d).*t(omega) + ...
            (0.5i*S(omega,lambda,tau).*t(omega))./q(omega));
    
    domg = omega_max / np;
    omg = transpose(domg : domg : 2 * omega_max);
    
    
    eomg = eb(omg, lambda, tau, z, d);
%     eomg(isnan(eomg)) = 0;
%     eomg(omg>2 * pi * 6e12) = 0.72 *  eomg(omg>2 * pi * 6e12);
    plot(omg/(2 * pi), abs(eomg).^2)
    ethz_omega = [zeros(4 * np, 1); eomg; zeros(2 * np, 1)];
    ethz_t = ifftshift(fft(fftshift(ethz_omega)));
    [~, max_ind ] =  max((real(ethz_t)));

    
    nf = length(ethz_omega);
    delta_t = 1/(nf * domg / (2 * pi));
    timeValues = transpose(1 : nf) * delta_t - delta_t * nf /2;

    time_ps = (timeValues - timeValues(max_ind))*1e12;
end


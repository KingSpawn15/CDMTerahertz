% 
% 

% 
% x = linspace(700, 800, 100);
% 
% y = materials.ngopt_znte(x * 1e-9);
% 
% nu = linspace(0,6,100);
% nTHz = materials.nTHz_znte(nu);
% 
% plot(nu,nTHz)

clear all


tau = 150e-15;
omega = 2 * pi * 1e12;
z = 2e-3;
d = 1e-3;
lambda = 800e-9;

% q0(omega, lambda)
% S(omega,lambda,tau)
% r(omega)
% t(omega)
% L(omega,lambda,z,d)
% 
% epl(omega, lambda, tau, z, d)
% emin(omega, lambda, tau, z, d)
% ef(omega, lambda, tau, z, d)
% eb(omega, lambda, tau, z, d)

%%





% domg = omega_max / np;
% omg = transpose(domg : domg : 2 * omega_max);
% 
% ethz_omega = [zeros(4 * np, 1); eb(omg, lambda, tau, z, d); zeros(2 * np, 1)];
% ethz_t = ifftshift(fft(fftshift(ethz_omega)));
% 
% nf = length(ethz_omega);
% 
% delta_t = 1/(nf * domg / (2 * pi));
% 
% timeValues = transpose(1 : nf) * delta_t - delta_t * nf /2;

%%
% [max_val, max_ind ] =  max(abs(real(ethz_t)));
close all;
omega_max = 2 * pi * 6e12;
np = 2001;

materials = opticalresponse;
nopt = @(lambda) materials.nopt_znte(1e6 * lambda);
nTHz = @(omega) materials.nTHz_znte(omega / (2 * pi) * 1e-12);
ngopt = @(lambda) materials.ngopt_znte(lambda * 1e6);


lambda_arr = [740, 770, 822, 860, 1560]*1e-9;
[time_ps1, ethz_t1] = electric_field_time(lambda_arr(1), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps2, ethz_t2] = electric_field_time(lambda_arr(2), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps3, ethz_t3] = electric_field_time(lambda_arr(3), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps4, ethz_t4] = electric_field_time(lambda_arr(4), tau, z, d, ngopt, nTHz, nopt, np, omega_max);
[time_ps5, ethz_t5] = electric_field_time(lambda_arr(5), tau, z, d, ngopt, nTHz, nopt, np, omega_max);

plot_scaling = 0.7;

plot(time_ps1, real(ethz_t1)/max(real(ethz_t3))*plot_scaling + 0.5 * 4 , ...
    time_ps2, real(ethz_t2)/max(real(ethz_t3))*plot_scaling+ 0.5 * 3,...
     time_ps3, real(ethz_t3)/max(real(ethz_t3))*plot_scaling+ 0.5 * 2,...
      time_ps4, real(ethz_t4)/max(real(ethz_t3))*plot_scaling+ 0.5 * 1,...
       time_ps5, real(ethz_t5)/max(real(ethz_t3))*plot_scaling+ 0.5 * 0,...
    'LineWidth',2);
xlim([-2,25])


% figure;
% 
% tau_arr = [300, 200, 150, 100, 50]*1e-15;
% [time_ps1, ethz_t1] = electric_field_time(lambda_arr(3), tau_arr(1), z, d, np, omega_max);
% [time_ps2, ethz_t2] = electric_field_time(lambda_arr(3), tau_arr(2), z, d, np, omega_max);
% [time_ps3, ethz_t3] = electric_field_time(lambda_arr(3), tau_arr(3), z, d, np, omega_max);
% [time_ps4, ethz_t4] = electric_field_time(lambda_arr(3), tau_arr(4), z, d, np, omega_max);
% [time_ps5, ethz_t5] = electric_field_time(lambda_arr(3), tau_arr(5), z, d, np, omega_max);
% 
% plot_scaling = 0.3;
% 
% plot(time_ps1, real(ethz_t1)/max(real(ethz_t3))*plot_scaling + 0.5 * 4 , ...
%     time_ps2, real(ethz_t2)/max(real(ethz_t3))*plot_scaling+ 0.5 * 3,...
%      time_ps3, real(ethz_t3)/max(real(ethz_t3))*plot_scaling+ 0.5 * 2,...
%       time_ps4, real(ethz_t4)/max(real(ethz_t3))*plot_scaling+ 0.5 * 1,...
%        time_ps5, real(ethz_t5)/max(real(ethz_t3))*plot_scaling+ 0.5 * 0,...
%     'LineWidth',2);
% xlim([-2,25])


%%


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
%     q0 = @(omega, lambda) omega .* materials.ngopt_znte(lambda * 1e6) / C;
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
    
    eb = @(omega, lambda, tau, z, d) emin(omega, lambda, tau, 0, d) .* t(omega) + ...
        (0.5i * S(omega, lambda, tau) .* t(omega)) ./ (exp(1i * z * qv(omega)) .* q(omega));

%     omega_max = 2 * pi * 2e12;
%     np = 1001;

    domg = omega_max / np;
    omg = transpose(domg : domg : 2 * omega_max);
    
    ethz_omega = [zeros(4 * np, 1); ef(omg, lambda, tau, z, d); zeros(2 * np, 1)];
    ethz_t = ifftshift(fft(fftshift(ethz_omega)));
    [~, max_ind ] =  max(abs(real(ethz_t)));

    nf = length(ethz_omega);
    delta_t = 1/(nf * domg / (2 * pi));
    timeValues = transpose(1 : nf) * delta_t - delta_t * nf /2;

    time_ps = (timeValues - timeValues(max_ind))*1e12;
end

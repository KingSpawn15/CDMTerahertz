function [time_ps, ethz_t, omg, eomg] = electric_field_time(lambda, tau, z, d, ngopt, nTHz, nopt, np, omega_max)


    C = 3e8;
    mu0 = 4 * pi * 1e-7;


    Ii = @(omega, tau) tau * sqrt(pi) * exp(-(tau^2 * omega.^2)/2);


    alpha = @(lambda)  (2 * pi / lambda) * imag(nopt(lambda));
    q = @(omega) omega .* nTHz(omega) / C;

%     q0 = @(omega, lambda) omega .* ngopt(lambda) / C  ;
    q0 = @(omega, lambda) real(omega .* ngopt(lambda) / C) + 1i * alpha(lambda) ;
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
    
   
     eb = @(omega, lambda, tau, z, d) exp(1i*z*qv(omega)).*(emin(omega,lambda,tau,0,d).*t(omega) + ...
            (0.5i*S(omega,lambda,tau).*t(omega))./q(omega));
    
    domg = omega_max / np;
    omg = transpose(domg : domg : 2 * omega_max);
    
    
    eomg = eb(omg, lambda, tau, z, d);
    plot(omg/(2 * pi), abs(eomg).^2)
    ethz_omega = [zeros(6 * np, 1); eomg; zeros(4 * np, 1)];
    ethz_t = ifftshift(fft(fftshift(ethz_omega)));
    [~, max_ind ] =  max((real(ethz_t)));

    
    nf = length(ethz_omega);
    delta_t = 1/(nf * domg / (2 * pi));
    timeValues = transpose(1 : nf) * delta_t - delta_t * nf /2;

    time_ps = (timeValues - timeValues(max_ind))*1e12;
end
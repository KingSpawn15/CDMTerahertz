function [TOR, ZOR, EOR] = electric_field_rectification(params_rectification, E_max_rectification, delay_or_pd_ps)
    % close all;
    % omega = 2 * pi * 1e12;
    
    tau = params_rectification.tau ;
    lambda = params_rectification.lambda ;
    d = params_rectification.d ;
    sigma_z = params_rectification.sigma_z ;
    z = params_rectification.z ;

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
    
    tt_t = time_ps - delay_or_pd_ps;


    ethz_t = real(ethz_t(tt_t < 10 & tt_t >-10));
    tt_t = tt_t(tt_t < 10 & tt_t >-10);
    vec_z = (-20 : 0.05 : 20 ).'*sigma_z;

    ethz_t = ethz_t/max(ethz_t);
%     E_max_rectification = (1.631e6) * 1.65;
    
    ethz_t = ethz_t .* E_max_rectification;
    
    [TOR, ZOR] = ndgrid(tt_t, vec_z);
    EOR = repmat(ethz_t, 1, length(vec_z)) .* exp(-ZOR.^2 / (2 * sigma_z^2));
    EOR = EOR.' ;
    ZOR = ZOR .* 1e6;
end
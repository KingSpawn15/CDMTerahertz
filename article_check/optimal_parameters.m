function fields = optimal_parameters()
    delete(gcp('nocreate'));
    parpool(6);
    
    % miscellaneous params
    e_w = linspace(-5,5,181);
    
    % photodember_parameters
    pump_power_nj = 10;
    laser_spot_size_fwhm = 40e-6;
    fitting_parameter_EPD = (1.26/6.34) * 1.2;
    
    % rectification parameters
    
    params_rectification.tau = 30e-15;
    params_rectification.lambda = 800e-9;
    params_rectification.d  = .5e-3;
    params_rectification.sigma_z = 50e-6;
    params_rectification.z = -1e-6;
    E_max_rectification = (1.631e6) * 1.65;
    delay_or_pd_ps = 0.1;
    
    % photodember field calculation
    eels_photodember = setup_parameters_eels_photodember(pump_power_nj, laser_spot_size_fwhm);
    [TPD, ZPD, EPD] = electric_field_photodember(eels_photodember, fitting_parameter_EPD);


    % Rectification field calculation
    
    [TOR, ZOR, EOR] = electric_field_rectification(params_rectification, E_max_rectification, delay_or_pd_ps);
    [T, Z, EPD, EOR] = interpolate_field(TOR, ZOR, EOR, TPD, ZPD, EPD);

    fields.EOR = EOR;
    fields.EPD = EPD;
    fields.e_w = e_w;
    fields.T = T;
    fields.Z = Z;
    fields.eels_obj = eels_photodember;
    
end

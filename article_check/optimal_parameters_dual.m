function fields = optimal_parameters_dual(weight, spot_size, shift, TPD_ni, ZPD_ni, EPD_ni, eels_photodember)
    

%     delete(gcp('nocreate'));
%     parpool(6)

%     pump_power_nj = 10;
%     laser_spot_size_fwhm = 40e-6;
%     fitting_parameter_EPD = (1.26/6.34) * 1.2;
%     eels_photodember = setup_parameters_eels_photodember(pump_power_nj, laser_spot_size_fwhm);
%     [TPD, ZPD, EPD] = electric_field_photodember(eels_photodember, fitting_parameter_EPD);
    % miscellaneous params
    e_w = linspace(-5,5,181);
    

   
    
    % rectification parameters
    
    params_rectification.tau = 30e-15;
    params_rectification.lambda = 800e-9;
    params_rectification.d  = .5e-3;
    params_rectification.sigma_z = 50e-6;
    params_rectification.z = -1e-6;
    E_max_rectification = (1.631e6) * 1.65;
    delay_or_pd_ps = shift;
    
%     % photodember field calculation


    % Rectification field calculation
    spot_size_fwhm = spot_size;
    [tmeep, zmeep, field_or_meep_0, field_or_meep_90] = rectification_meep_dual(spot_size_fwhm);
    
    tmeep = tmeep - delay_or_pd_ps;
    [TOR, ZOR] = ndgrid(tmeep, zmeep);
    EOR_0 = field_or_meep_0.';
    EOR_90 = field_or_meep_90.';
    
    mEOR = max(EOR_0(:));
    EOR_0  = EOR_0 / mEOR * E_max_rectification * weight;
    EOR_90 = EOR_90 / mEOR * E_max_rectification * weight;
    %%
%     [TOR, ZOR, EOR] = electric_field_rectification(params_rectification, E_max_rectification, delay_or_pd_ps);
    [~, ~, ~, EOR_0] = interpolate_field(TOR, ZOR, EOR_0, TPD_ni, ZPD_ni, EPD_ni);
    [T, Z, EPD, EOR_90] = interpolate_field(TOR, ZOR, EOR_90, TPD_ni, ZPD_ni, EPD_ni);

    fields.EOR_0 = EOR_0;
    fields.EOR_90 = EOR_90;
    fields.EPD = EPD;
    fields.e_w = e_w;
    fields.T = T;
    fields.Z = Z;
    fields.eels_obj = eels_photodember;
    
end

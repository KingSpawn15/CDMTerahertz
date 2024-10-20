function optimization_power_weights(weights, EOR_10, EPD, T, Z)
    
    max_exp = [0.3500    0.4000    0.5500    0.8000    0.8000    1.1000    1.8000    2.5000    3.3000    3.6500    4.2000];
%     weights = [0.1541    0.2360    0.4436    0.6615    0.5648    0.5694    0.7773    0.9870    1.1000    1.4222    1.5135];
    
    for ii = 1:11
        disp(ii)
        [t0_vec,eels_comb{ii}] = utils_spectrum.calculate_spectrum_from_fields_coarse(EOR_10 * weights(ii) + EPD{ii}, T, Z * 1e-6);
    end

    max_sim = [];
    for ii = 1:11
        max_sim = [max_sim,max(abs(eels_comb{ii}))];
    end
    
    
end

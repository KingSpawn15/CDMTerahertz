function [tc, xc, field_laser_profile_x, field_laser_profile_y, field_laser_profile_z] = get_fields_rectification(spot_size_fwhm)
%     load('meep_results\saved_matrices_meep\triple_test_ez_3\field_ez50.0_fsshift0.3_ps.mat')
%     load('meep_results\saved_matrices_meep\rectification_eps_12\field_ez50_fsshift0.1_ps.mat')
    load('meep_results\saved_matrices_meep\rectification_eps_12\field_ez_with_derivative_50_fsshift0.1_ps.mat')

    
    xc =  - zstep * size(e_or_x,2) / 2 : zstep:  zstep * size(e_or_x,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_or_x,1);
    field_laser_profile_x = convolve_laser_profile(e_or_x, xc, spot_size_fwhm).';
    field_laser_profile_y = convolve_laser_profile(e_or_y, xc, spot_size_fwhm).';
    field_laser_profile_z = convolve_laser_profile(e_or_z, xc, spot_size_fwhm).';

end

%%
function field_laser_profile = convolve_laser_profile(field, xc, spot_size_fwhm)
    
    sigma = spot_size_fwhm / sqrt(8 * log(2));
    intensity = exp(-(xc).^2/(2 * sigma^2)) ;

    for i = 1:size(field, 1)
        field_laser_profile(i, :) = conv_line(field(i, :), intensity, xc);
    end

end

function conv_fi = conv_line(field, intensity, xc)

        padded_intensity = [zeros(1,length(xc)), intensity,zeros(1,length(xc))];
        
        padded_field = [zeros(1,length(xc)), field,zeros(1,length(xc))];
        conv_fi = conv(padded_intensity, padded_field, 'same');
        conv_fi = conv_fi(length(xc)+1:2*length(xc));
end

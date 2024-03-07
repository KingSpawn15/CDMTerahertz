function [tc, xc, field_laser_profile_0, field_laser_profile_90] = rectification_meep_dual(spot_size_fwhm)
    load('article_check/saved_matrices_meep/set6/field_ez50.0_fs.mat')
    
    %%
%     close all
%     imagesc(e_or_0.')
    
    %%
  
    xc =  - zstep * size(e_or_90,2) / 2 : zstep:  zstep * size(e_or_90,2) / 2 - zstep;
    tc = tstep : tstep : tstep * size(e_or_90,1);
    field_laser_profile_0 = convolve_laser_profile(e_or_0, xc, spot_size_fwhm);
    field_laser_profile_90 = convolve_laser_profile(e_or_90, xc, spot_size_fwhm);
%     
%     close all
%     imagesc(tc, xc, field_laser_profile.')

    
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

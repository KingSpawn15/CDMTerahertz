function [interact_v] = eels_pdor(eels,method)
    
switch method

    case 'photodember'

        loss_spectrum_parameters.method = 'photodember';
        loss_spectrum_parameters.interaction_gain_factor_rectification = 0;
        loss_spectrum_parameters.interaction_gain_factor_photodember = 1;
        interact_v = eels.interaction_v(loss_spectrum_parameters);

    case 'rectification'
        loss_spectrum_parameters.method = 'rectification';
        loss_spectrum_parameters.interaction_gain_factor_rectification = 1;
        loss_spectrum_parameters.interaction_gain_factor_photodember = 0;
        interact_v = eels.interaction_v(loss_spectrum_parameters);
    otherwise
        disp('method rect/pd only')
end

end


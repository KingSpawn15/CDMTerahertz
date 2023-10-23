function [TPD, ZPD, EPD] = electric_field_photodember(eels_photodember, fitting_parameter_EPD)
    
    interact_v_pd = utils.eels_pdor(eels_photodember,'photodember');
    [TPD, ZPD] = ndgrid(eels_photodember.discretization.t * 1e12,eels_photodember.discretization.z * 1e6);
    EPD = utils_spectrum.derivative_f_dz(interact_v_pd, eels_photodember.discretization.z, ...
        eels_photodember.discretization.t) * fitting_parameter_EPD;
    EPD(isnan(EPD)) = 0;

end
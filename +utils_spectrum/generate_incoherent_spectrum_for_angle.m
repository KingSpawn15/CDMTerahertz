function [e_w, t_w, psi_incoherent] = generate_incoherent_spectrum_for_angle(EPD, EOR, T, Z, theta, eels_photodember, e_w)
    
    [t_w, eels_spectra] = utils_spectrum.calculate_spectrum_from_fields(EPD -cos(2 * theta * pi / 180) * EOR, T, Z*1e-6);
    psi_assemb = utils_spectrum.spectrum_to_coherent_eels_mod(e_w, t_w, eels_spectra);
    w = eels_photodember.electron.incoherent_gaussian_blurring_window(e_w,t_w);
    psi_incoherent =  eels_photodember.incoherent_convolution(psi_assemb, w, t_w, e_w);

end
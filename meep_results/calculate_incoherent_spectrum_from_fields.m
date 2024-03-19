function [t_w , psi_incoherent] = calculate_incoherent_spectrum_from_fields(EZ, T, Z, e_w)

    [t_w, eels_spectra] = utils_spectrum.calculate_spectrum_from_fields(EZ, T, Z*1e-6);
    psi_assemb = utils_spectrum.spectrum_to_coherent_eels_mod(e_w, t_w, eels_spectra * 1e2 * 3);
    eels_photodember = setup_parameters_eels_photodember(10, 40e-6);
    w = eels_photodember.electron.incoherent_gaussian_blurring_window(e_w,t_w);
    psi_incoherent =  eels_photodember.incoherent_convolution(psi_assemb, w, t_w, e_w);

end
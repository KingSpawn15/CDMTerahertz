function [t0_vec, eels] = calculate_spectrum_from_fields(EOR, TOR, ZOR)

EOR = EOR.';
max_z = 50e-6;
velec = 0.7 * 3 * 10^(8-12);

electron_z = (-5 : 0.01: 5).' .* max_z;
electron_t0 = electron_z / velec;
t0_vec = (-3:0.02:3).';
eels = zeros(length(t0_vec),1);


length_t0_vec = length(t0_vec);

parfor ind = 1: length_t0_vec
  
    et_l = interp2(TOR', ZOR', EOR', electron_t0 + t0_vec(ind), electron_z, "cubic",0);
    eels(ind) = trapz(electron_z,et_l);

end


end
% function [t0_vec, eels] = calculate_spectrum_from_fields(EOR, TOR, ZOR)
% 
% EOR = EOR.';
% max_z = 50e-6;
% velec = 0.7 * 3 * 10^(8-12);
% 
% electron_z = (-5 : 0.1: 5).' .* max_z;
% electron_t0 = electron_z / velec;
% t0_vec = (-3:0.02:3).';
% eels = zeros(length(t0_vec),1);
% 
% 
% length_t0_vec = length(t0_vec);
% 
% parfor ind = 1: length_t0_vec
%   
%     et_l = interp2(TOR', ZOR', EOR', electron_t0 + t0_vec(ind), electron_z, "linear",0);
%     eels(ind) = trapz(electron_z,et_l);
% 
% end
% 
% 
% end

function [t0_vec, eels] = calculate_spectrum_from_fields_coarse(EOR, TOR, ZOR)

EOR = EOR.';
max_z = 50e-6;
velec = 0.7 * 3 * 10^(8-12);

electron_z = (-3 : 0.01: 3).' .* max_z;
electron_t0 = electron_z / velec;
t0_vec = (-1.5:0.05:1.5).';
eels = zeros(length(t0_vec),1);


length_t0_vec = length(t0_vec);

for ind = 1: length_t0_vec
  
    et_l = interp2(TOR', ZOR', EOR', electron_t0 + t0_vec(ind), electron_z, "linear",0);
%     eels(ind) = trapz(electron_z,et_l);
    eels(ind) = trapz(et_l) * 0.01 * max_z;

end


end
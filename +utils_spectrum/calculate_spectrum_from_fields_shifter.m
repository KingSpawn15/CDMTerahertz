function [t0_vec, eels] = calculate_spectrum_from_fields_shifter(EOR, TOR, ZOR)



EOR = EOR.';
% max_z = 90e-6;
velec = 0.7 * 3 * 10^(8-12);
nz = size(EOR,2);

dz = ZOR(1,2) - ZOR(1,1);
dt = TOR(2,1) - TOR(1,1);

pad_EOR = [EOR*0;EOR;EOR*0];

for i = 1:nz
    
    pad_EOR(:,i) = circshift(pad_EOR(:,i),-round((i-nz/2) * dz / (velec)/dt));
end

eels_padded = sum(pad_EOR,2);
t_points = -TOR(end,1)+dt:dt:2*TOR(end,1);
%     pad_EOR(:,i) = circshift(pad_EOR(:,i),[round(-nz/2 + i * dz / (2 * velec)/dt),0]);

% size(EOR)
% electron_z = (-3 : 0.01: 3).' .* max_z;
% electron_t0 = electron_z / velec;
t0_vec = (-2.3:0.01:2.3).';
% eels = zeros(length(t0_vec),1);
eels = interp1(t_points, eels_padded, t0_vec) * dz;


% length_t0_vec = length(t0_vec);
% 
% parfor ind = 1: length_t0_vec
%   
%     et_l = interp2(TOR', ZOR', EOR', electron_t0 + t0_vec(ind), electron_z, "linear",0);
%     eels(ind) = trapz(electron_z,et_l);
% 
% end


end
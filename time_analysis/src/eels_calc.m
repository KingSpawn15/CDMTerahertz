function [T, Z, ET, t0_vec, eels] = eels_calc(ethz, t, sigma_z, velec)

vec_z = (-10 : 0.05 : 10 ).'*sigma_z;
[T, Z] = ndgrid(t, vec_z);
ET = repmat(ethz, 1, length(vec_z)) .* exp(-Z.^2 / (2 * sigma_z^2));

t0_vec = (-2:0.02:3);
eels = [];
for t0 = t0_vec

    zz_l = (-5 : 0.1: 5).' .* sigma_z;
    tt_l = zz_l / velec - t0;
    et_l = interp2(T', Z', ET', tt_l, zz_l);
%     et_l(isnan(et_l)) = 0;
    eels = [eels;trapz(zz_l,et_l)];
    
end


t0_vec = t0_vec.';




end


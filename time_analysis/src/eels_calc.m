function [T, Z, ET, t0_vec, eels] = eels_calc(ethz, t, sigma_z, velec)

vec_z = (-20 : 0.05 : 20 ).'*sigma_z;
[T, Z] = ndgrid(t, vec_z);
ET = repmat(ethz, 1, length(vec_z)) .* exp(-Z.^2 / (2 * sigma_z^2));
% ET(abs(Z)>100e-6) =0;
t0_vec = (-5:0.02:5);
eels = zeros(length(t0_vec),1);
ind = 1;
for t0 = t0_vec

    zz_l = (-5 : 0.1: 5).' .* sigma_z;
    tt_l = zz_l / velec + t0;
    et_l = interp2(T', Z', ET', tt_l, zz_l);
%     et_l(isnan(et_l)) = 0;
    eels(ind) = trapz(zz_l,et_l);
    ind = ind + 1;
end


t0_vec = t0_vec.';




end


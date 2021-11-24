function [fwhm_t , fwhm_e] = calculate_marginals(w, e_w, t_w)
%CALCULATE_MARGINALS Summary of this function goes here
%   Detailed explanation goes here


proj_t = trapz(e_w,w,2);
proj_e = trapz(t_w,w,1);

[~,max_ind] = max(proj_t);
e1 = find(proj_t >= 0.5*proj_t(max_ind),1,'first');
e2 = find(proj_t >= 0.5*proj_t(max_ind),1,'last');
fwhm_t = t_w(e2) - t_w(e1);

[~,max_ind] = max(proj_e);
e1 = find(proj_e >= 0.5*proj_e(max_ind),1,'first');
e2 = find(proj_e >= 0.5*proj_e(max_ind),1,'last');
fwhm_e = e_w(e2) - e_w(e1);



end


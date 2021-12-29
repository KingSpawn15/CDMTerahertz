function err = error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory)
%ERROR_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
params.interaction_gain_factor_photodember = factor_pd_theory;
params.interaction_gain_factor_rectification = factor_or_theory;
params.theta_pol_degree = angle_theory;
[psi_incoherent, e_w, t_w] = optimization.wrapper_polarization(params);

angle = angle_measurement;
[psi_measurement, e_w_measurement,...
    t_w_measurement] = measurement_plots.data_measurement(angle);

window_e_min = -2.5;
window_e_max = 2.5;

window_t_min = -1;

window_t_max = 1;




[E_W , T_W] = meshgrid(e_w, t_w);
[E_W_measurement , T_W_measurement] = meshgrid(e_w_measurement, t_w_measurement);

window = (heaviside(E_W - window_e_min) - heaviside(E_W - window_e_max)).* ...
    (heaviside(T_W - window_t_min) - heaviside(T_W - window_t_max));

psi_theory_sub = interp2(E_W, T_W, window .* psi_incoherent, E_W_measurement, T_W_measurement,...
    'nearest',0);

err = psi_measurement - psi_theory_sub;


end


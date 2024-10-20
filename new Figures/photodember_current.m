pump_power_nj = 10;
laser_spot_size_fwhm = 40e-6;
eels_photodember = setup_parameters_eels_photodember(pump_power_nj, laser_spot_size_fwhm);

nexc = eels_photodember.laser.excited_carriers(eels_photodember.material.alpha, eels_photodember.material.hnew);
omega_y_0 = eels_photodember.material.calculate_omega_y_0(nexc);

t_ps = -0.5:0.01:3;
sigma_t = 50e-15 / (2*sqrt(2 * log(2)));
t_sec = t_ps * 1e-12;

% gaussian_pulse = exp(-t_sec .^2 / (2 * sigma_t^2));



current_pd = sin(omega_y_0 * t_ps * 1e-12) .* exp(-eels_photodember.material.gamma * t_ps * 1e-12 / 2);
current_pd(t_ps < 0) = 0;

gaussian_pulse_ps = exp(-t_ps .^2 / (2 * (sigma_t*1e12)^2)) ;
gaussian_pulse_ps = gaussian_pulse_ps / max(gaussian_pulse_ps) * max(current_pd);
%%
close all
h = plot(t_ps, current_pd, t_ps, gaussian_pulse_ps);
xlim([-.3,1.5]);
axis square
set(gca, 'FontSize', 14);
set(h, {'LineWidth'}, {1;0.5});
exportgraphics(gcf, ['new Figures/results/', 'section_jPD.png'],'resolution', 300);

%%
function y = convolve_same_length(x, h)
    % Ensure the vectors are of the same size
    if length(x) ~= length(h)
        error('Vectors must be of the same length');
    end
    
    % Zero padding only at the end
    x_padded = [x, zeros(1, length(x)-1)];
    h_padded = [h, zeros(1, length(h)-1)];
    
    % Perform convolution
    y_full = conv(x_padded, h_padded);
    
    % Return the result of the same length as the input vectors
    y = y_full(1 + 52:length(x) + 52);
end


function current_pd_time = current_pd_convolved(pulse_fwhm_fs,t_ps, current_pd)
    sigma_t = pulse_fwhm_fs * 1e-15 / (2*sqrt(2 * log(2)));
    t_sec = t_ps * 1e-12;
    
    gaussian_pulse = exp(-t_sec .^2 / (2 * sigma_t^2));
    
    current_pd_time = convolve_same_length(gaussian_pulse, current_pd);
    current_pd_time = current_pd_time / max(current_pd_time) * max(current_pd);
end


function [alpha_interp] = interpolate_absorption(data, lambda_nm)
    % Extract data from the struct
    lambda_data = data.lambda_nm;
%     n_data = data.n;
    k_data = data.k;

    % Perform linear interpolation
%     n_interp = interp1(lambda_data, n_data, lambda_nm, 'linear', 'extrap');
    k_interp = interp1(lambda_data, k_data, lambda_nm, 'linear', 'extrap');
    
    Cnm = 3e17;
    f = Cnm ./ lambda_nm;
    omega = 2 * pi * f;
    alpha_interp = 2 .* omega ./ Cnm .* k_interp;
end
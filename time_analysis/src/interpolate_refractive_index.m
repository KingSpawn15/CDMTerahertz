function ngopt = interpolate_refractive_index(data, lambda_nm)
    % Extract data from the struct
    lambda_data = data.lambda_nm;
    n_data = data.n;
    k_data = data.k;

    % Perform linear interpolation
    n_interp = interp1(lambda_data, n_data, lambda_nm, 'linear', 'extrap');
    k_interp = interp1(lambda_data, k_data, lambda_nm, 'linear', 'extrap');
    ngopt = n_interp + 1i * k_interp;
end
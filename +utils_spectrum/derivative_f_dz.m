function dfdz = derivative_f_dz(f_values, z_values, t_values)
    % Check if the input matrices have compatible sizes
    if size(f_values, 1) ~= length(z_values) || size(f_values, 2) ~= length(t_values)
        error('Input matrix dimensions do not match with vector sizes.');
    end

    % Calculate the derivative using the finite difference method
    dz = z_values(2) - z_values(1); % Assuming uniform grid spacing
    dfdz = diff(f_values, 1, 1) ./ dz;

    % Pad the dfdz matrix with NaNs to maintain the same size as f_values
    dfdz = [dfdz; NaN(1, size(dfdz, 2))];
end

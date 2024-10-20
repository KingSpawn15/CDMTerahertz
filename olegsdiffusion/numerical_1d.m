main()

function main
    % Simulation parameters
    L = 2;               % Length of the spatial domain
    T = 3;               % Total simulation time
    Nx = 100;            % Number of spatial grid points
    Nt = 100;            % Number of time steps
    mu = 0.1;            % Electron mobility
    D = 0.01;            % Diffusion constant
    tr = 1;              % Recombination time
    q = 1;         % Charge of electron
    epsilon_0 = 1; % Permittivity of free space

    % Define the source function
    S_func = @(x, t) exp(-7 * x) .* exp(-((t - 1).^2) / (2 * 0.05^2));  % Example source function
    
    % Call the function to solve the density equations
    [n_e, n_h, E, x_grid, t_grid] = solve_density_eqns_reflect(L, T, Nx, Nt, mu, D, tr, q, epsilon_0, S_func);
end

function [n_e, n_h, E, x_grid, t_grid] = solve_density_eqns_reflect(L, T, Nx, Nt, mu, D, tr, q, epsilon_0, S_func)
    % Discretize the spatial and time domains
    dx = L / (Nx - 1);
    dt = T / (Nt - 1);
    x_grid = linspace(0, L, Nx);  % Domain [0, L]
    t_grid = linspace(0, T, Nt);
    
    % Initialize the electron and hole densities
    n_e = zeros(Nx, Nt);
    n_h = zeros(Nx, Nt);
    
    % Set initial conditions for electron and hole densities based on the source term at t = 0
    S_initial = S_func(x_grid, 0);    % Source term at t = 0
    n_e(:, 1) = S_initial(:);         % Initialize electron density from the source term
    n_h(:, 1) = S_initial(:);         % Initialize hole density from the source term
    
    % Finite difference matrices
    D2x = (diag(-2 * ones(Nx, 1)) + diag(ones(Nx-1, 1), 1) + diag(ones(Nx-1, 1), -1)) / dx^2;

    % Apply Neumann BC (zero-flux) at the boundaries (for simplicity)
    D2x(1, 2) = 2 / dx^2;       % Reflecting BC at x = 0 (Neumann)
    D2x(end, end-1) = 2 / dx^2; % Reflecting BC at x = L (if needed)
    
    % Initialize the electric field
    E = zeros(Nx, Nt);
    
    % Time-stepping loop
    for n = 1:Nt-1
        t = t_grid(n);
        
        % Compute source term at the current time step
        Sx_t = S_func(x_grid, t);
        Sx_t = Sx_t(:);  % Ensure Sx_t is a column vector of size Nx

        % Calculate the electric field E
        E(:, n) = compute_electric_field(n_e(:, n), n_h(:, n), q, epsilon_0, dx);
        
        % Update n_e with the electric field
        gradE_n_e = mu * E(:, n);  % Use electric field directly for update
        n_e(:, n+1) = n_e(:, n) + dt * (gradE_n_e + D * D2x * n_e(:, n) - (n_e(:, n) .* n_h(:, n)) / tr + Sx_t);
        
        % Reflecting boundary at x=0: enforce n_e(1) = n_e(2)
        n_e(1, n+1) = n_e(2, n+1);
        
        % Hole density update with the influence of the electric field
        n_h(:, n+1) = n_h(:, n) + dt * (- (n_e(:, n) .* n_h(:, n)) / tr + Sx_t);
        
        % Reflecting boundary at x=0: enforce n_h(1) = n_h(2)
        n_h(1, n+1) = n_h(2, n+1);
    end
    
    % Plot the result (optional)
    figure;
    subplot(3, 1, 1);
    imagesc(t_grid, x_grid, n_e);
    colorbar;
    xlabel('Time');
    ylabel('Position (x)');
    title('Electron Density n_e(x,t)');
    
    subplot(3, 1, 2);
    imagesc(t_grid, x_grid, n_h);
    colorbar;
    xlabel('Time');
    ylabel('Position (x)');
    title('Hole Density n_h(x,t)');
    
    subplot(3, 1, 3);
    imagesc(t_grid, x_grid, E);
    colorbar;
    xlabel('Time');
    ylabel('Position (x)');
    title('Electric Field E(x,t)');
end

function E = compute_electric_field(n_e, n_h, q, epsilon_0, dx)
    % Calculate the electric field E based on the Poisson equation
    rho = q * (n_e - n_h);
    E = zeros(size(n_e));  % Initialize electric field
    
    % Using central difference to solve for E
    E(2:end-1) = E(2:end-1) + (rho(2:end-1) / epsilon_0) * (dx^2);
    
    % Reflecting boundary conditions for electric field
    E(1) = E(2);      % Reflecting boundary at x = 0
    E(end) = E(end-1); % Reflecting boundary at x = L (if needed)
    E = E*0;
end

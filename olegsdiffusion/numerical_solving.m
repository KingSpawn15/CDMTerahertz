clear all
close all
% Parameters
nc = 4.6e3;          % number density n0
alpha = 7;          % Decay constant
gamma = 3.3;        % Decay constant for source term
sigma_y = 20;
D = 0.1;              % Diffusion coefficient
t_r = 100;          % Recombination time (you need to define this)
Lx = 5;             % Length in x-direction
Ly = 200;           % Length in y-direction
T = 3;              % Total time
Nx = 49;            % Number of grid points in x-direction (reduced)
Ny = 500;           % Number of grid points in y-direction
Nt = 1000;          % Number of time steps (increased for more stability)
t0 = 0.1;
sigma_t = 0.02;
mu = 4;

dx = Lx / (Nx - 1); % Grid spacing in x-direction
dy = Ly / (Ny - 1); % Grid spacing in y-direction
dt = T / Nt;        % Time step size

% Stability condition check
if D * dt / dx^2 > 0.5 || D * dt / dy^2 > 0.5
    error('The solution is unstable. Reduce dt or increase dx and dy.');
end

% Create grid using meshgrid
x = linspace(-Lx, Lx, Nx);
y = linspace(-Ly, Ly, Ny);
[X, Y] = meshgrid(x, y);  % Create a meshgrid for x and y
% Initial condition function
eps = ones(size(X));
eps(X<0) = eps(X<0)*12;
eps = eps';
n_e = nc * exp(-alpha * X) .* exp(-(Y.^2) / (2 * sigma_y^2)) .* exp(-(t0)^2 / (2 * sigma_t^2));  % Initial condition at t = 0
n_e(X<0) = 0;
n_e = n_e';                % Transpose for correct orientation
n_h = n_e;                 % Initialize hole density from electron density
div_neE  = n_e * 0;
phi_old = div_neE + 1e-6;

matEy = [];
%%
% Create a figure for visualization
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);  % [left, bottom, width, height]

% Time-stepping loop
for n = 1:Nt

    n_e_new = n_e; % Create a new matrix for the updated electron densities
    n_h_new = n_h; % Create a new matrix for the updated hole densities

    % Calculate source term
    source_term = nc * (1 / sqrt( 2 * pi * sigma_t^2)) * exp(-alpha * X) .* exp(-(Y.^2) / (2 * sigma_y^2)) .* exp(-(n * dt-t0)^2 / (2 * sigma_t^2));
    source_term(X<0) = 0;
    source_term = source_term';
    
    % Calculate second-order derivatives using finite difference (vectorized)
    d2_n_e_dx2 = (n_e(3:Nx, 2:Ny-1) - 2 * n_e(2:Nx-1, 2:Ny-1) + n_e(1:Nx-2, 2:Ny-1)) / dx^2;
    d2_n_e_dy2 = (n_e(2:Nx-1, 3:Ny) - 2 * n_e(2:Nx-1, 2:Ny-1) + n_e(2:Nx-1, 1:Ny-2)) / dy^2;
    
    % Update electron density n_e (without the loop, using matrix operations)
    dn_dt = zeros(size(n_e_new));
    dn_dt(2:Nx-1, 2:Ny-1) = (D * (d2_n_e_dx2 + d2_n_e_dy2) + source_term(2:Nx-1, 2:Ny-1) - ...
         (n_h(2:Nx-1, 2:Ny-1) .* n_e(2:Nx-1, 2:Ny-1) / t_r) + mu * div_neE(2:Nx-1, 2:Ny-1) + ...
         - n_e(2:Nx-1, 2:Ny-1) * gamma);
    n_e_new(2:Nx-1, 2:Ny-1) = n_e(2:Nx-1, 2:Ny-1) + dt * dn_dt(2:Nx-1, 2:Ny-1);
    
    n_h_new(2:Nx-1, 2:Ny-1) = n_h(2:Nx-1, 2:Ny-1) + dt * ...
        (-n_h(2:Nx-1, 2:Ny-1) .* n_e(2:Nx-1, 2:Ny-1) / t_r + ...
        + source_term(2:Nx-1, 2:Ny-1));
    
    % Reflecting boundary condition at x = 0 (assuming x corresponds to the index 25)
    n_e_new((Nx+1)/2, :) = n_e_new((Nx+1)/2 + 1, :);  
    n_h_new((Nx+1)/2, :) = n_h_new((Nx+1)/2 + 1, :);    
    
    % Absorbing boundary conditions at y = -Ly (lower y-bound) and y = Ly (upper y-bound)
    n_e_new(1 : (Nx-1)/2, :) = 0;    % Absorbing at both extremes of x (y = -Ly, y = Ly)
    n_e_new([1, Nx], :) = 0;  % Absorbing condition at top and bottom boundaries in x
    n_e_new(:, [1, Ny]) = 0;    % Absorbing at both extremes of y

    n_h_new(1 : (Nx-1)/2, :) = 0;    % Apply the same for n_h
    n_h_new(:, [1, Ny]) = 0;    % Apply the same for n_h
    n_h_new([1, Nx], :) = 0;  % Absorbing condition at top and bottom boundaries in x

    % Update n_e and n_h for the next time step
    n_e = n_e_new;
    n_h = n_h_new;

    % compute divergence for next iteration
    [Ex, Ey, div_neE, phi] = compute_divergence_and_field(n_e, n_h, x, y, phi_old, eps);
    [jx, jy] = compute_current(n_e, Ex, Ey, D, mu, dx, dy);
    phi_old = phi;
    
    matEy = horzcat(matEy,Ey(20,:)');

    if mod(n,20) == 0
        params.n = n;
        params.dt = dt;
        params.x = x;
        params.y = y;
        params.n_e = n_e;
        params.n_h = n_h;
        params.source_term = source_term;
        params.Ex = Ex;
        params.Ey = Ey;
        params.X = X;
        params.Y = Y;
        params.nc = nc;
        params.jx = jx;
        params.jy = jy;
        params.dn_dt = dn_dt;
        plotResults(params);
    end



end



%%
function [j_x, j_y] = compute_current(n_e, E_x, E_y, D, mu, dx, dy)
    % Compute the gradients of electron density
    q = 1.6e-19;
    [d_n_e_dy, d_n_e_dx] = gradient(n_e, dy, dx);
    
    % Calculate the current densities
    j_x = -q * D * d_n_e_dx + q * mu * n_e .* E_x;  % Current density in the x-direction
    j_y = -q * D * d_n_e_dy + q * mu * n_e .* E_y;  % Current density in the y-direction
end

function [Ex, Ey, div_neE, phi] = compute_divergence_and_field(n_e, n_h, x, y, phi, eps)
    % Constants
    epsilon_0 = 8.854e-18; % Permittivity of free space
    

    % Compute electric field from densities
    tol = 1e-5;
    max_iter = 5000;
    [Ex, Ey, phi, converged] = poisson_solver(n_e, n_h, x, y, epsilon_0, phi, tol, max_iter, eps);
    
    if converged == false
        warning("poisson solver not converged");
    end

    dx = x(2) - x(1);
    dy = y(2) - y(1);
    [d_n_e_dy, d_n_e_dx] = gradient(n_e, dy, dx);
    
    % Compute divergence of the electric field E
    [dE_x_dy, dE_x_dx] = gradient(Ex, dy, dx);
    [dE_y_dy, dE_y_dx] = gradient(Ey, dy, dx);

    % Apply the product rule to compute divergence
    div_neE = Ex .* d_n_e_dx + Ey .* d_n_e_dy +  n_e .* (dE_x_dx + dE_y_dy);
end

function [E_x, E_y, phi, converged] = poisson_solver(n_e, n_h, x, y, epsilon_0, phi, tol, max_iter, eps)
    % Define the grid
%     [X, Y] = meshgrid(x, y);
    
    % Define the charge density
    qe = -1.60217662e-19; % Elementary charge in Coulombs
    rho = qe * (n_e - n_h);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    % Initialize the potential with the initial guess
    V = phi;
    
    % Set boundary conditions (example: V = 0 at the boundaries)
    V(:, 1) = 0; % Left boundary
    V(:, end) = 0; % Right boundary
    V(1, :) = 0; % Bottom boundary
    V(end, :) = 0; % Top boundary
    
    % Iterative solver (Gauss-Seidel method)
    converged = false;
    for iter = 1:max_iter
        V_old = V;
        
        for i = 2:size(n_h,1)-1
            for j = 2:size(n_h,2)-1
                % Average permittivity at half-grid points
                epsilon_x1 = (eps(i,j) + eps(i+1,j)) / 2;
                epsilon_x2 = (eps(i,j) + eps(i-1,j)) / 2;
                epsilon_y1 = (eps(i,j) + eps(i,j+1)) / 2;
                epsilon_y2 = (eps(i,j) + eps(i,j-1)) / 2;
                
                % Update phi at point (i, j) using neighboring values
                V(i,j) = (1 / (epsilon_x1/dx^2 + epsilon_x2/dx^2 + epsilon_y1/dy^2 + epsilon_y2/dy^2)) * ...
                    (epsilon_x1 * V(i+1,j) / dx^2 + ...
                     epsilon_x2 * V(i-1,j) / dx^2 + ...
                     epsilon_y1 * V(i,j+1) / dy^2 + ...
                     epsilon_y2 * V(i,j-1) / dy^2 + ...
                     rho(i,j) / epsilon_0);
            end
        end

        if max(max(abs(V - V_old))) < tol
            converged = true;
            disp(iter)
            break;
        end
    end
    
    % Calculate the electric field components
    [E_y, E_x] = gradient(-V, dy, dx);
    
    
    % Return the potential and convergence flag
    phi = V;
end

function plotResults(params)
    % Plotting every 20 time steps
    % Create a subplot for n_e


    subplot(2, 3, 1); 
    imagesc(params.x, params.y, params.n_e'); 
    xlabel('x');
    ylabel('y');
    zlabel('n_e(x, y, t)');
    title(['Electron Density at t = ', num2str(params.n * params.dt)]);
    colorbar;
    hold on;  % To add the white line on top of the plot
    line([0 0], ylim, 'Color', 'w', 'LineWidth', 0.5,'LineStyle','-'); % Draw vertical white line at x=0
    hold off;

    % Create a subplot for n_h
    subplot(2, 3, 2); 
    imagesc(params.x, params.y, params.n_h'); 
    xlabel('x');
    ylabel('y');
    zlabel('n_h(x, y, t)');
    title(['Hole Density at t = ', num2str(params.n * params.dt)]);
    colorbar;
    hold on;  % To add the white line on top of the plot
    line([0 0], ylim, 'Color', 'w', 'LineWidth', 0.5,'LineStyle','-'); % Draw vertical white line at x=0
    hold off;

    % Create a subplot for source_term
    subplot(2, 3, 3); 
    imagesc(params.x, params.y, params.source_term'); 
    xlabel('x');
    ylabel('y');
    zlabel('Source Density (x, y, t)');
    title(['Source Density at t = ', num2str(params.n * params.dt)]);
    clim([0, params.nc]);
    colorbar;
    hold on;  % To add the white line on top of the plot
    line([0 0], ylim, 'Color', 'w', 'LineWidth', 0.5,'LineStyle','-'); % Draw vertical white line at x=0
    hold off;

    % Create a subplot for Ex
    subplot(2, 3, 4); 
    imagesc(params.x, params.y, params.jx'); 
    xlabel('x');
    ylabel('y');
    zlabel('Jx');
    title(['Jx at t = ', num2str(params.n * params.dt)]);
    colorbar;
    hold on;  % To add the white line on top of the plot
    line([0 0], ylim, 'Color', 'w', 'LineWidth', 0.5,'LineStyle','-'); % Draw vertical white line at x=0
    hold off;

    % Create a subplot for Ex
    subplot(2, 3, 5); 
    imagesc(params.x, params.y, params.jy'); 
    xlabel('x');
    ylabel('y');
    zlabel('Jy');
    title(['Jy at t = ', num2str(params.n * params.dt)]);
    colorbar;
    hold on;  % To add the white line on top of the plot
    line([0 0], ylim, 'Color', 'w', 'LineWidth', 0.5,'LineStyle','-'); % Draw vertical white line at x=0
    hold off;


    % Create a subplot for Ey
%     subplot(2, 3, 6); 
%     surf(params.X, params.Y, params.dn_dt', 'EdgeColor', 'none'); 
%     xlabel('x');
%     ylabel('y');
%     zlabel('Ey');
%     title(['Ey at t = ', num2str(params.n * params.dt)]);
%     view(2);  % Set the view to top-down
%     colorbar;
%     hold on;  % To add the white line on top of the plot
%     line([0 0], ylim, 'Color', 'w', 'LineWidth', 0.5,'LineStyle','-'); % Draw vertical white line at x=0
%     hold off;

    drawnow;
end


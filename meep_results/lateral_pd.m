% Your parameters
N0 = 1;
gamma = 3.3;
v = 7.6e-1;
sigma = 40/(sqrt(8 * log (2)));
l = v / gamma;

% Define the range of x and t
x = linspace(-100, 100, 500); % 1000 points between -100 and 100
t = linspace(-5, 10, 1000); % 1000 points between -5 and 10

% Create a grid of x and t values
[X, T] = meshgrid(x, t);

% Heaviside function in MATLAB is 'heaviside'
H = double(T >= 0);

% Your expression
Z = (exp(-(X.^2./(2*(l.*T.*v + sigma^2)))) .* N0 .* X .* H) ./ (sqrt(2*pi) .* sqrt((l.*T.*v)./sigma + sigma) .* (l.*T.*v + sigma^2));

% Define your t0 and sigma_t for the Gaussian function
t0 = 1; % Your defined t0
sigma_t = 50e-5; % Your defined sigma_t

% Gaussian function
gaussian_func = exp(-(t - t0).^2 / (2*sigma_t^2));

% Pad the Gaussian function with zeros
pad_size = 0;
gaussian_func_padded = [zeros(1, pad_size), gaussian_func, zeros(1, pad_size)];

% Flip the Gaussian function for convolution
% gaussian_func_padded = flip(gaussian_func_padded);

% Convolute each column of Z with the padded Gaussian function
Z_conv = zeros(size(Z)); % Initialize the convoluted matrix
for i = 1:size(Z, 2)
    conv_res = conv(Z(:, i), gaussian_func_padded, 'same');
    % Trim the convolution result back to the original size
    Z_conv(:, i) = conv_res(2*pad_size+1:length(Z(:, i))+2*pad_size);
end

% Flip the convoluted matrix back to the original time direction
Z_conv = flip(Z_conv, 1);

% Create the imagesc plot with t on the x-axis and x on the y-axis
close all
imagesc(t, x, Z_conv');
colorbar; % Add a colorbar to the plot
xlabel('t');
ylabel('x');
title('Imagesc plot after convolution');

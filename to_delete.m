Fs = 44100;         % Sampling frequency
T = 1/Fs;           % Sampling period
t = -0.5:T:0.5;     % Time vector
L = length(t);      % Signal length

X =sin(2 * pi * 2000 * t);
n = 2^nextpow2(L);
Y = (fft(X,n));

f = Fs*(0:(n/2))/n;
P = abs(Y/n).^2;

plot(f,P(1:n/2+1)) 
title("Gaussian Pulse in Frequency Domain")
xlabel("f (Hz)")
ylabel("|P(f)|^2")

% 
% X = 1/(0.4*sqrt(2*pi))*(exp(-t.^2/(2*(0.1*1e-3)^2)));
% % Plot the pulse in the time domain.
% 
% plot(t,X)
% title("Gaussian Pulse in Time Domain")
% xlabel("Time (t)")
% ylabel("X(t)")
% axis([-1e-3 1e-3 0 1.1]) 
% 
% n = 2^nextpow2(L);
% % Convert the Gaussian pulse to the frequency domain.
% 
% Y = fftshift(fft(X,n));
% % Define the frequency domain and plot the unique frequencies.
% 
% f = Fs*(0:(n/2))/n;
% P = abs(Y/n).^2;
% 
% plot(f,P(1:n/2+1)) 
% title("Gaussian Pulse in Frequency Domain")
% xlabel("f (Hz)")
% ylabel("|P(f)|^2")


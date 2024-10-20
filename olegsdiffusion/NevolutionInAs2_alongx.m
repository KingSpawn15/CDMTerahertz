% electron expansion modeling with x coordinate and exp(-alpha * x) initial condition

clear all, close all,

dt = 0.01; % in ps
toffst = 0.2;
pulseD = 0.1;
t = 0.00:dt:3; % time vector
[unkn, ttot] = size(t);

xlim = 2;  % Limit for x in microns
xcel = 0.05; % Spatial resolution in microns
x = 0:xcel:xlim; % space vector for x coordinate
[unkn, xcells] = size(x);

alpha = 7;  % Decay constant for exp(-alpha * x)
TotEl = 1000; % Total number of electrons generated, NOTE: 1D model
niniF = exp(-alpha * x);  % Exponential initial condition in x
niniTempF = exp((-(t-toffst).^2)/pulseD^2);  % Time-dependent factor
ampEl = TotEl/(sum(niniF)*sum(niniTempF)); 
nini = niniF * ampEl;

NxaEvol(1:xcells, 1:ttot) = 0;
NxaEvolPump(1:xcells, 1:ttot) = 0;

for iii = 1:20
    NxaEvolPump(1:xcells, iii) = nini * niniTempF(iii);
end

Cur(1:xcells, 1:ttot) = 0;

% Assuming velocity is constant, 10 microns/ps
elecVumpps = 10;
vel = elecVumpps / xcel * dt;  % in cells per time interval

figure(1),
subplot(1,2,1)
plot(x, NxaEvol(1:xcells, 1));
subplot(1,2,2)
plot(t(1:100), niniTempF(1,1:100));

% for tstep = 2:ttot
%     ElDistTemp(1) = 0;
%     ElDistTemp(2:(xcells+1)) = NxaEvol(1:xcells, (tstep-1));
%     ElDistTemp((xcells+2)) = 0;
%     
%     NxaEvol(1:xcells, tstep) = NxaEvolPump(1:xcells, tstep); 
%     NxaEvol(1:xcells, tstep) = NxaEvol(1:xcells, tstep) + NxaEvol(1:xcells, (tstep-1)) * (1 - 2 * vel);  % initial number minus number leaving the cell
%     NxaEvol(1:xcells, tstep) = NxaEvol(1:xcells, tstep) + vel * (ElDistTemp(1:(xcells)) + ElDistTemp(3:(xcells+2)))';  % number entering the cell
%     Cur(1:xcells, tstep) = vel * (-ElDistTemp(1:(xcells)) + ElDistTemp(3:(xcells+2)))';
%     CurDer(1:xcells, tstep) = Cur(1:xcells, tstep) - Cur(1:xcells, (tstep-1));
% end

for tstep = 2:ttot
    % Temporary variable to store the electron distribution
    ElDistTemp(1) = NxaEvol(1, (tstep-1));  % Reflecting boundary at x = 0
    ElDistTemp(2:(xcells+1)) = NxaEvol(1:xcells, (tstep-1));  % Shifted array for handling boundaries
    ElDistTemp((xcells+2)) = 0;  % Initialize for the boundary at the other end (absorbing boundary)

    % Update the electron density for the current time step
    NxaEvol(1:xcells, tstep) = NxaEvolPump(1:xcells, tstep);  % Pumping term
    NxaEvol(1:xcells, tstep) = NxaEvol(1:xcells, tstep) + NxaEvol(1:xcells, (tstep-1)) * (1 - 2 * vel);  % Number minus those leaving the cell
    NxaEvol(1:xcells, tstep) = NxaEvol(1:xcells, tstep) + vel * (ElDistTemp(1:(xcells)) + ElDistTemp(3:(xcells+2)))';  % Adding those entering the cell

    % Reflecting boundary condition at x = 0:
    NxaEvol(1, tstep) = NxaEvol(2, tstep);  % Reflecting the density at x = 0 back into the system

    % Absorbing boundary condition at the other end (x = 5, i.e., xcells):
    NxaEvol(xcells, tstep) = 0;  % Absorb all particles at the boundary (set electron density to 0)

    % Current calculation
    Cur(1:xcells, tstep) = vel * (-ElDistTemp(1:(xcells)) + ElDistTemp(3:(xcells+2)))';
    CurDer(1:xcells, tstep) = Cur(1:xcells, tstep) - Cur(1:xcells, (tstep-1));

    % Ensure the current at x = 0 reflects as well
    Cur(1, tstep) = -Cur(2, tstep);  % Reflect the current at the boundary

    % Absorbing boundary condition for current at x = 5 (xcells):
    Cur(xcells, tstep) = 0;  % Ensure no current is leaving at the absorbing boundary
end



figure(2),

subplot(1,3,1)
plot(x, NxaEvolPump(1:xcells, 1), '-k'),
hold on,
plot(x, NxaEvolPump(1:xcells, 6), '-b'),
hold on,
plot(x, NxaEvolPump(1:xcells, 10), '-r'),
hold on,
plot(x, NxaEvolPump(1:xcells, 15), '-g'),
hold off;

subplot(1,3,2)
plot(x, NxaEvol(1:xcells, 1), '-k'),
hold on,
plot(x, NxaEvol(1:xcells, 100), '-b'),
hold on,
plot(x, NxaEvol(1:xcells, 200), '-r'),
hold off;

subplot(1,3,3)
plot(x, Cur(1:xcells, 1), '-k'),
hold on,
plot(x, Cur(1:xcells, 100), '-b'),
hold on,
plot(x, Cur(1:xcells, 200), '-r'),
hold off;

tsh=0:dt:1;
[TT XX] = meshgrid(tsh, x);

figure(3),

subplot(2,2,1)
surf(TT, XX, NxaEvolPump(1:xcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 0 xlim]);
xlabel('time (ps)'),
ylabel('x (um)'),
title('Intensity')

subplot(2,2,2)
surf(TT, XX, NxaEvol(1:xcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 0 xlim]);
xlabel('time (ps)'),
ylabel('x (um)'),
title('Elec. Density')

subplot(2,2,3)
surf(TT, XX, Cur(1:xcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 0 xlim]);
xlabel('time (ps)'),
ylabel('x (um)'),
title('Current')

subplot(2,2,4)
surf(TT, XX, CurDer(1:xcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 0 xlim]);
xlabel('time (ps)'),
ylabel('x (um)'),
title('dj/dt')

figure(4),

subplot(1,4,1)
imagesc(NxaEvolPump(1:xcells, 1:100));

subplot(1,4,2)
imagesc(NxaEvol(1:xcells, 1:100));

subplot(1,4,3)
imagesc(Cur(1:xcells, 1:100));

subplot(1,4,4)
imagesc(CurDer(1:xcells, 1:100));

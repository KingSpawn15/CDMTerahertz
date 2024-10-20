% electon expansion modelling

clear all, close all,

dt = 0.01 % in ps
toffst = 0.1;
pulseD = 0.025;
t = 0.00:dt:3; % time vector
[unkn, ttot] = size(t);

zalim = 50;
zcel = 1;
za = (-1*zalim):zcel:(zalim); % space vector
[unkn, zcells] = size(za);

wid = 20; % Gaussian width of the optical beam in microns
TotEl = 1000; % total number of electrons generated NOTE: 1D model
niniF = exp((-za.^2)/wid^2);
niniTempF = exp((-(t-toffst).^2)/pulseD^2);
ampEl = TotEl/(sum(niniF)*sum(niniTempF)); 
nini = niniF*ampEl;

NzaEvol(1:zcells, 1:ttot) = 0;

NzaEvolPump(1:zcells, 1:ttot) = 0;

for iii = 1:20;
    NzaEvolPump(1:zcells, iii) = nini*niniTempF(iii);
end;

Cur(1:zcells, 1:ttot) = 0;

% assuming velocity is constant and it is 10 u/ps (the initial velocity can
% be higher, can be calculated from the kinetic energy given the photon
elecVumpps = 10;
vel = elecVumpps/zcel * dt; % in cells per time interval
% then multipying the vel by Npercell, one can find the number of electrons
% leaving the cell in one time step;

figure(1), 
subplot(1,2,1)
plot(za, NzaEvol(1:zcells, 1));
subplot(1,2,2)
plot(t(1:100), niniTempF(1,1:100));


for tstep = 2:ttot,
    ElDistTemp(1) = 0;
    ElDistTemp(2:(zcells+1)) = NzaEvol(1:zcells, (tstep-1));
    ElDistTemp((zcells+2)) = 0;
    
    NzaEvol(1:zcells, tstep) = NzaEvolPump(1:zcells, tstep); 
    NzaEvol(1:zcells, tstep) = NzaEvol(1:zcells, tstep) + NzaEvol(1:zcells, (tstep-1)) * (1 - 2 * vel); %initial number minus the number leaving the cell
    NzaEvol(1:zcells, tstep) = NzaEvol(1:zcells, tstep) + vel*(ElDistTemp(1:(zcells))+ElDistTemp(3:(zcells+2)))'; %prev line plus the number entering the cell
    Cur(1:zcells, tstep) = vel*(-ElDistTemp(1:(zcells))+ElDistTemp(3:(zcells+2)))';
    CurDer(1:zcells, tstep) = Cur(1:zcells, tstep) - Cur(1:zcells, (tstep-1));  
end;

figure(2),

subplot(1,3,1)
plot(za, NzaEvolPump(1:zcells, 1), '-k'),
hold on,
plot(za, NzaEvolPump(1:zcells, 6), '-b'),
hold on,
plot(za, NzaEvolPump(1:zcells, 10), '-r'),
hold on,
plot(za, NzaEvolPump(1:zcells, 15), '-g'),
hold off;

subplot(1,3,2)
plot(za, NzaEvol(1:zcells, 1), '-k'),
hold on,
plot(za, NzaEvol(1:zcells, 100), '-b'),
hold on,
plot(za, NzaEvol(1:zcells, 200), '-r'),
hold off;

subplot(1,3,3)
plot(za, Cur(1:zcells, 1), '-k'),
hold on,
plot(za, Cur(1:zcells, 100), '-b'),
hold on,
plot(za, Cur(1:zcells, 200), '-r'),
hold off;

tsh=0:dt:1;
[TT ZZ] = meshgrid(tsh, za);

figure(3),

subplot(2,2,1)
surf(TT, ZZ, NzaEvolPump(1:zcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 -50 50]);
xlabel('time (ps)'),
ylabel('z (um)'),
title('Intensity')


subplot(2,2,2)
surf(TT, ZZ, NzaEvol(1:zcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 -50 50]);
xlabel('time (ps)'),
ylabel('z (um)'),
title('Elec. Density')


subplot(2,2,3)
surf(TT, ZZ, Cur(1:zcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 -50 50]);
xlabel('time (ps)'),
ylabel('z (um)'),
title('Current')


subplot(2,2,4)
surf(TT, ZZ, CurDer(1:zcells, 1:101));
shading("interp");
view(2);
axis([0 0.6 -50 50]);
xlabel('time (ps)'),
ylabel('z (um)'),
title('dj/dt')



figure(4),

subplot(1,4,1)
imagesc(NzaEvolPump(1:zcells, 1:100));

subplot(1,4,2)
imagesc(NzaEvol(1:zcells, 1:100));

subplot(1,4,3)
imagesc(Cur(1:zcells, 1:100));

subplot(1,4,4)
imagesc(CurDer(1:zcells, 1:100));


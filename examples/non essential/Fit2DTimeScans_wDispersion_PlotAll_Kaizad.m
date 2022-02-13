clc; close all; clear all;

%%
InputArgs = [1.4e-02,7.1e-02,1.3e+13,2.7e-01,-8.1e-01];

PulseEnergyGainFactor = InputArgs(1);
PhiGainFactor = InputArgs(2);
Gamma = InputArgs(3);
GammaFactor = InputArgs(4);
Phase = InputArgs(5);

%%
% Fitted parameters
ElectronTotalEnergy = 1.1;
ElectronTotalTime = 3.6e-13;
LaserSpotFWHM = 40e-6;

ElectronTimeCoherentFWHM = 50e-15;
Theta = -7*pi/180;

%% Physical constants
c = 299792458;%[m/s];
e = 1.60217662e-19;%[C]
h = 6.62607004e-34;%[J*s]
hbar = h/(2*pi);

%% System parameters
%Coherent parameters
ElectronTimeCoherentSigma = ElectronTimeCoherentFWHM./(2*sqrt(2*log(2)));%[s]

%Incoherent parameters
ElectronTimeIncoherentSigma = ElectronTotalTime/(2*sqrt(2*log(2)));%[s]    
ElectronEnergyIncoherentSigma = ElectronTotalEnergy/(2*sqrt(2*log(2)));%[eV]

%
v = 0.7*c;%[m/s]
C1 = (-1./(1i*hbar*v))*e;

%% Time variable
Fs = 2.4e15;%[Hz]
dt = 1/Fs;%[s]
L = 2.4e4;
domega = Fs/L*(2*pi);%[Hz]
t = (-L/2:L/2).*dt;%[s]
omega = (-L/2:L/2)*domega;%[Hz]
energy = omega.*hbar/e;%[eV]

%% Time delay
ddt = 10e-15;
deltat = -2.5e-12:ddt:2.5e-12;

%% Incoherent broadening weighting function W(E,t)
SubSampleFactor = 60;
eW = energy(1:SubSampleFactor:end);%[eV]
tW = deltat*1e12;%[ps]

SigmaX = ElectronTimeIncoherentSigma*1e12;
SigmaY = ElectronEnergyIncoherentSigma;

a = cos(Theta)^2/(2*SigmaX^2) + sin(Theta)^2/(2*SigmaY^2);
b = (sin(2*Theta)/4)*((1/SigmaY^2)-(1/SigmaX^2));
c = sin(Theta)^2/(2*SigmaX^2) + cos(Theta)^2/(2*SigmaY^2);

[X,Y] = meshgrid(tW,eW);
W = exp(-(a*X.^2 + 2*b*X.*Y + c*Y.^2));
W = W';

figure(10)
imagesc(eW,tW,W)
xlabel('Energy [eV]')
ylabel('Time [ps]')
axis square
drawnow

ProjT = trapz(eW,W,2);
ProjE = trapz(tW,W,1);

[~,MaxInd] = max(ProjT);
e1 = find(ProjT >= 0.5*ProjT(MaxInd),1,'first');
e2 = find(ProjT >= 0.5*ProjT(MaxInd),1,'last');
FWHM_T = tW(e2) - tW(e1);

[~,MaxInd] = max(ProjE);
e1 = find(ProjE >= 0.5*ProjE(MaxInd),1,'first');
e2 = find(ProjE >= 0.5*ProjE(MaxInd),1,'last');
FWHM_E = eW(e2) - eW(e1);

fprintf('Projected energy FWHM: %.3f [eV]\nProjected time FWHM: %.3f [ps]\n',FWHM_E,FWHM_T);

%% z variable
z = (-1:2e-2:1)*1e-4;%[m]
deltaz = z(2)-z(1);

[Omega,Z] = meshgrid(omega,z);

%% Scan parameters
PulseEnergy = PulseEnergyGainFactor.*10.000.*1e-9;%[J]

x0 = 0;
y0 = -1e-6;

%% Calculate electric potential

Phi = CalcElectricPotential_wXprimeZprimeInt_wRetPotential_womega_xy(PulseEnergy,t,z,x0,y0,Gamma,GammaFactor,LaserSpotFWHM,Phase);
Phi = Phi.*PhiGainFactor;



Phi = movmean(movmean(Phi,3,1),10,2);

figure(20)
imagesc(t,z,Phi)
xlabel('t [s]')
ylabel('z [m]')
drawnow

%% Time domain Fourier transform
PhiFFT = fftshift(fft(Phi,length(t),2),2);
PhiFFT = PhiFFT./max(omega);

%% beta(omega)
PhiZ = PhiFFT.*exp(-1i*Omega.*Z/v);
beta = sum(PhiZ,1).*deltaz;

%% f(t)
parfor TimeInd = 1:length(t)
    I(TimeInd) = C1.*trapz(omega,2.*real(exp(1i.*omega.*t(TimeInd)).*beta));    
end

%% Sweep time delay - Psi(E,deltat)
Deltat = repmat(deltat',1,length(t));
T = repmat(t,length(deltat),1);
Ishift = circshift(I,ceil(length(t)/2));
ft = exp(-Ishift);
FT = repmat(ft,length(deltat),1);
PsiCoherent = FT.*exp(-(T-Deltat).^2/(2*ElectronTimeCoherentSigma.^2));

PsiCoherent = fftshift(fft(PsiCoherent,length(t),2),2);

PsiCoherent = (abs(PsiCoherent)).^2;
PsiCoherent = PsiCoherent./trapz(energy,PsiCoherent,2);

params_m.PulseEnergy = PulseEnergy;
params_m.Gamma = Gamma;
params_m.GammaFactor = GammaFactor;
params_m.LaserSpotFWHM = LaserSpotFWHM;
params_m.Phase = Phase;
params_m.Phi = Phi;
params_m.PsiCoherent = PsiCoherent;

save('params_m.mat','params_m');

figure(30)
imagesc(energy,deltat.*1e12,PsiCoherent)
title('|\psi_{Coherent}(E)|^2','FontWeight','Normal')
xlabel('Energy [eV]')
ylabel('\Deltat [ps]')
colorbar
colormap jet
axis square
drawnow

%% Sub-sample the model time scan and normalize
TimeScan_Model = PsiCoherent(:,1:SubSampleFactor:end);
TimeScan_Model = TimeScan_Model./trapz(eW,TimeScan_Model,2);

%% Incoherent convolution
TimeScan_Model_Sum = zeros(size(TimeScan_Model));

WCutOff = 0.01*max(W(:));

parfor tInd = 1:length(tW)
    for eInd = 1:length(eW)

        if W(tInd,eInd) < WCutOff
            continue
        end

        TimeScan_Model_Sum = TimeScan_Model_Sum + W(tInd,eInd).*...
            circshift(circshift(TimeScan_Model,-ceil(length(tW)/2)+tInd,1),-ceil(length(eW)/2)+eInd,2);

    end
end

TimeScan_Model_HR = TimeScan_Model_Sum;
TimeScan_Model_HR = TimeScan_Model_HR./trapz(eW,TimeScan_Model_HR,2);
TimeScan_Model_HR = TimeScan_Model_HR./max(TimeScan_Model_HR(:));

%% Plot
figure(40)
imagesc(eW,tW,TimeScan_Model_HR)
title('\psi_{Model}','FontWeight','Normal')
xlabel('Energy [eV]')
ylabel('\Deltat [ps]')
xlim([-3 3])
ylim([-0.5 1.2])
colorbar
colormap jet
axis square
drawnow
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 18;
ax.XTick = -3:1.5:3;
ax.YTick = -0.5:0.5:1;

phi_incoherent_m = TimeScan_Model_HR;
save('phi_incoherent_m','phi_incoherent_m');

%%
clc; close all; clear all;

debug = true;

%% Input arguments
if debug
    InputArgs = [0.8,150,20,-5*pi/180,1e-1,1];
end

ElectronTotalEnergy = InputArgs(1);
ElectronTotalTime = InputArgs(2)*1e-15;
ElectronTimeCoherentFWHM = InputArgs(3)*1e-15;
Theta = InputArgs(4);
PhiGainFactor = InputArgs(5);
PulseEnergyGainFactor = InputArgs(6);

fprintf('Input arguments:\n\t 1. ElectronTotalEnergy: %.3f [eV]\n\t 2. ElectronTotalTime: %.3f [fs]\n\t 3. ElectronTimeCoherentFWHM: %.3f [fs]\n\t 4. Theta: %.3f [deg]\n\t 5. PhiGainFactor: %.3f\n\t 6. PulseEnergyGainFactor: %.3f\n',InputArgs.*[1 1 1 180/pi 1 1]);

%% Predetermined parameters
Theta_Pol = 90.*(pi/180);

PulseEnergy = PulseEnergyGainFactor.*1e-10;%[J]

x0 = 0;
y0 = -1e-6;

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

v = 0.7*c;% Electron velocity [m/s]
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

%% z variable
z = (-1:1e-2:1)*1e-4;%[m]
deltaz = z(2)-z(1);

%% Incoherent broadening weighting function W(E,t)
SubSampleFactor = 60;
eW = energy(1:SubSampleFactor:end);
tW = deltat*1e12;%[ps]

SigmaT = ElectronTimeIncoherentSigma*1e12;
SigmaE = ElectronEnergyIncoherentSigma;

a = cos(Theta)^2/(2*SigmaT^2) + sin(Theta)^2/(2*SigmaE^2);
b = (sin(2*Theta)/4)*((1/SigmaE^2)-(1/SigmaT^2));
c = sin(Theta)^2/(2*SigmaT^2) + cos(Theta)^2/(2*SigmaE^2);

[TW,EW] = meshgrid(tW,eW);
W = exp(-(a*TW.^2 + 2*b*TW.*EW + c*EW.^2));
W = W';

figure(1)
imagesc(eW,tW,W)
xlabel('Energy [eV]')
ylabel('Time [ps]')
axis square
drawnow

%% Calculate the marginals of W
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

%% Calculate electric potential
% Here you can choose which potential to calculate (or use both and sum)
% Phi_OR - Optical Rectificatio
% Phi_PD - Photo-Dember

Phi_OR = CalcElectricPotential_OpticalRectification_wRetPotential(PulseEnergy,t,z,x0,y0,Theta_Pol,v);
% Phi_PD = CalcElectricPotential_wXprimeZprimeInt_wRetPotential(PulseEnergy,t,z,x0,y0);

%%
Phi = Phi_OR.*PhiGainFactor;

Phi = movmean(movmean(Phi,3,1),10,2);

figure(2)
imagesc(t,z,Phi)
xlabel('t [s]')
ylabel('z [m]')
drawnow

%% Time domain Fourier transform
PhiFFT = fftshift(fft(Phi,length(t),2),2);
PhiFFT = PhiFFT./max(omega);

%% beta(omega)
[Omega,Z] = meshgrid(omega,z);

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

%% Sub-sample the model time scan and normalize
TimeScan_Model = PsiCoherent(:,1:SubSampleFactor:end);
TimeScan_Model = TimeScan_Model./trapz(eW,TimeScan_Model,2);

figure(3)
imagesc(eW,tW,TimeScan_Model)
xlabel('Energy [eV]')
ylabel('\Deltat [ps]')
colorbar
colormap jet
axis square
drawnow

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

%%
TimeScan_Model = TimeScan_Model_Sum;
TimeScan_Model = TimeScan_Model./trapz(eW,TimeScan_Model,2);

%% Plot final model output TimeScan(E,deltat)
figure(4)
imagesc(eW,tW,TimeScan_Model)
xlabel('Energy [eV]')
ylabel('\Deltat [ps]')
colorbar
colormap jet
axis square
drawnow
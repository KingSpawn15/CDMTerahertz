% function [RSSError] = Fit2DTimeScans_OR_wDispersion(InputArgs)

ExperimentInd = 12;

% PulseEnergyGainFactor_OR = InputArgs(1);
% PhiGainFactor_OR = InputArgs(2);

PulseEnergyGainFactor_OR = 1;       % fitting parameter for pulse energy
PhiGainFactor_OR = 1e-2;            % fitting paramer potential for V

% Fitted parameters
ElectronTotalEnergy = 1.2;          % eV (FWHM)
ElectronTotalTime = 3.600e-13;      % pulse duration

ElectronTimeCoherentFWHM = 12e-15;  % incoherent broadening
Theta = -20/180*pi;                 % Angle of tilt in phase space

%% Physical constants
c = 299792458;%[m/s];
e = 1.60217662e-19;%[C]
h = 6.62607004e-34;%[J*s]
hbar = h/(2*pi);
hbar_eVs = h/e;

%% System parameters
%Coherent parameters
ElectronTimeCoherentSigma = ElectronTimeCoherentFWHM./(2*sqrt(2*log(2)));%[s]
ElectronEnergyCoherentSigma = (hbar_eVs/2)./ElectronTimeCoherentSigma;%[eV]
ElectronEnergyCoherentFWHM = ElectronEnergyCoherentSigma.*(2*sqrt(2*log(2)));

%Incoherent parameters
ElectronTimeIncoherentFWHM_P = sqrt(ElectronTotalTime^2 - (ElectronTimeCoherentFWHM./sqrt(2)).^2);
ElectronTimeIncoherentSigma_P = ElectronTimeIncoherentFWHM_P/(2*sqrt(2*log(2)));%[s]    
ElectronEnergyIncoherentFWHM_P = sqrt(ElectronTotalEnergy^2 - (ElectronEnergyCoherentFWHM./sqrt(2)).^2);
ElectronEnergyIncoherentSigma_P = ElectronEnergyIncoherentFWHM_P/(2*sqrt(2*log(2)));%[eV]

% 
v = 0.7*c;%[m/s]            % 200 keV
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
ddt = 5e-15;
deltat = -5e-12:ddt:5e-12;          % Time delay vector

%% z variable
z = (-1:1e-2:1)*1e-4;%[m]
deltaz = z(2)-z(1);

[Omega,Z] = meshgrid(omega,z);

%% Scan parameters

pulse_energy_experiment = 10e-9;
PulseEnergy = PulseEnergyGainFactor_OR*pulse_energy_experiment;
%PulseEnergy = PulseEnergy(ExperimentInd);

Theta_Pol = (-2:2:2).*(pi/180);

Length = length(Theta_Pol);

% Where the electron is passing. 
x0 = 0;
y0 = -1e-6;

for PolInd = 1:Length

    %% Calculate electric potential

    [Phi] = CalcElectricPotential_OpticalRectification_wRetPotential(PulseEnergy,t,z,x0,y0,...
        Theta_Pol(PolInd), v);
    % Phi is V
    Phi = Phi.*PhiGainFactor_OR;

    PhiAll(PolInd,:,:) = Phi;

    Phi = movmean(movmean(Phi,3,1),10,2);       % make smoother
    
    figure(10)
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
    FT = repmat(ft,length(deltat),1);  % time and time delay matrix
    PsiCoherent = FT.*exp(-(T-Deltat).^2/(2*ElectronTimeCoherentSigma.^2));

    PsiCoherent = fftshift(fft(PsiCoherent,length(t),2),2);

    PsiCoherent = (abs(PsiCoherent)).^2;
    PsiCoherent = PsiCoherent./trapz(energy,PsiCoherent,2);

    PsiCoherentAll(PolInd,:,:) = PsiCoherent;
    
    figure(20 + (PolInd-1))
    imagesc(energy,deltat.*1e12,PsiCoherent)
    title('|\psi_{Coherent}(E)|^2','FontWeight','Normal')
    xlabel('Energy [eV]')
    ylabel('\Deltat [ps]')
    colorbar
    colormap hot
    axis square
    drawnow
    
end

% end
clc; close all; clear all;

InputArgs = [1.4e-02,7.1e-02,1.3e+13,2.7e-01,-8.1e-01];
ExpInd = 1:12;

%%
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
% PulseEnergy = PulseEnergyGainFactor.*[0.017,0.030,0.054,0.100,0.170,0.300,0.540,1.000,1.770,3.000,5.400,10.000].*1e-9;%[J]
% PulseEnergy = PulseEnergyGainFactor.*[0.01:0.02:0.1, 0.2:0.2:1, 1.5:0.5:10].*1e-9;%[J]
PulseEnergy = PulseEnergyGainFactor.*[10].*1e-9;%[J]
Length = length(PulseEnergy);

x0 = 0;
y0 = -1e-6;

for PulseEnergyInd = 1:Length%ExpInd

    fprintf('Pulse energy No. %i\n',PulseEnergyInd); 
    
    pulseenergy = PulseEnergy(PulseEnergyInd);

    %% Calculate electric potential
    Phi = CalcElectricPotential_wXprimeZprimeInt_wRetPotential_womega_xy(pulseenergy,t,z,x0,y0,Gamma,GammaFactor,LaserSpotFWHM,Phase);
    Phi = Phi.*PhiGainFactor;
    
    % Case Phi is imaginary
    if (~isreal(Phi))
        RSSError = 1e10;
        return;
    end
    
    Phi = movmean(movmean(Phi,3,1),10,2);

    figure(200 + PulseEnergyInd)
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

    figure(300 + PulseEnergyInd)
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
    
    %%
%     load('TimeScansAll_Experiment.mat');
%     load('CrossSectionDeltatAll.mat');
% 
%     TimeScan_Experiment = TimeScansAll_Experiment{PulseEnergyInd};
%     TimeScan_Experiment = TimeScan_Experiment./trapz(Energy_Experiment,TimeScan_Experiment,2);
% 
%     for tInd = 1:length(TimeDelay_Experiment)
%         TimeDelayInd(tInd) = find(deltat*1e12 >= TimeDelay_Experiment(tInd),1,'first');
%     end
%     for eInd = 1:length(Energy_Experiment)
%         EnergyInd(eInd) = find(eW >= Energy_Experiment(eInd),1,'first');
%     end
%     dt_Experiment = TimeDelay_Experiment(2) - TimeDelay_Experiment(1);
%     dE_Experiment = Energy_Experiment(2) - Energy_Experiment(1);
    
    %%
%     TimeScan_Model = TimeScan_Model_Sum(TimeDelayInd,EnergyInd);
%     TimeScan_Model = TimeScan_Model./trapz(Energy_Experiment,TimeScan_Model,2);
    
    TimeScan_Model_HR = TimeScan_Model_Sum;
    TimeScan_Model_HR = TimeScan_Model_HR./trapz(eW,TimeScan_Model_HR,2);

    %%
%     TimeScan_Model = TimeScan_Model./max(TimeScan_Model(:));
%     TimeScan_Experiment = TimeScan_Experiment./max(TimeScan_Experiment(:));

    TimeScan_Model_HR = TimeScan_Model_HR./max(TimeScan_Model_HR(:));
    
    %% <Energy>
%     S_Model = trapz(Energy_Experiment,Energy_Experiment.*TimeScan_Model,2)./trapz(Energy_Experiment,TimeScan_Model,2);
%     S_Experiment = trapz(Energy_Experiment,Energy_Experiment.*TimeScan_Experiment,2)./trapz(Energy_Experiment,TimeScan_Experiment,2);
    
    S_Model_HR = trapz(eW,eW.*TimeScan_Model_HR,2)./trapz(eW,TimeScan_Model_HR,2);
    
    %% Align the model time scan to the experiment
%     MaxTimeShift = 20;
%     EnergyShift = -0.2:0.01:0.2;
%     deW = eW(2) - eW(1);
% 
%     for tInd = -MaxTimeShift:1:MaxTimeShift
%         for eInd = 1:length(EnergyShift)
%             
%             S_Model_Shift = circshift(S_Model,tInd) + EnergyShift(eInd);
%             Error(MaxTimeShift+1+tInd,eInd) = sqrt(sum( (abs(S_Model_Shift - S_Experiment)).^2 ));
%         
%         end
%     end
% 
%     [~,ErrorMinInd] = min(Error(:));
%     [ErrorMinIndT,ErrorMinIndE] = ind2sub(size(Error),ErrorMinInd);
%         
%     ErrorMinIndTAll(PulseEnergyInd) = ErrorMinIndT;
%     ErrorMinIndEAll(PulseEnergyInd) = ErrorMinIndE;
%     
%     S_Model = circshift(S_Model,ErrorMinIndT-(MaxTimeShift+1)) + EnergyShift(ErrorMinIndE);
%     TimeScan_Model = circshift(circshift(TimeScan_Model,ErrorMinIndT-(MaxTimeShift+1),1),round(EnergyShift(ErrorMinIndE)/dE_Experiment),2);
%     
%     S_Model_HR = circshift(S_Model_HR,round((ErrorMinIndT-(MaxTimeShift+1)).*(dt_Experiment/(ddt*1e12)))) + EnergyShift(ErrorMinIndE);
%     TimeScan_Model_HR = circshift(circshift(TimeScan_Model_HR,round((ErrorMinIndT-(MaxTimeShift+1)).*(dt_Experiment/(ddt*1e12))),1)...
%         ,round(EnergyShift(ErrorMinIndE)/deW),2);

    %%
%     S_Experiment_All(PulseEnergyInd,:) = S_Experiment;
%     TimeScan_Experiment_All(PulseEnergyInd,:,:) = TimeScan_Experiment;

%     S_Model_All(PulseEnergyInd,:) = S_Model;
%     TimeScan_Model_All(PulseEnergyInd,:,:) = TimeScan_Model;

    S_Model_HR_All(PulseEnergyInd,:) = S_Model_HR;
    TimeScan_Model_HR_All(PulseEnergyInd,:,:) = TimeScan_Model_HR;
    
    %% Plot
    figure(400 + PulseEnergyInd)
    subplot(1,2,1)
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
% 
%     subplot(1,2,2)
%     imagesc(Energy_Experiment,TimeDelay_Experiment,TimeScan_Experiment)
%     title('\psi_{Exp}','FontWeight','Normal')
%     xlabel('Energy [eV]')
%     ylabel('\Deltat [ps]')
%     xlim([-3 3])
%     ylim([-0.5 1.2])
%     colorbar
%     colormap jet
%     axis square
%     drawnow
%     box on
%     ax = gca;
%     ax.LineWidth = 2;
%     ax.FontSize = 18;
%     ax.XTick = -3:1.5:3;
%     ax.YTick = -0.5:0.5:1;
% 
%     %% Plot
%     figure(500 + PulseEnergyInd)
%     hold on
%     plot(tW,S_Model_HR,'LineWidth',2)
%     plot(TimeDelay_Experiment,S_Experiment,'LineWidth',2)
%     hold off
%     legend({'<E>_{Model}','<E>_{Exp}'})
%     xlabel('\Deltat [ps]')
%     xlim([-0.5 1.2])
%     axis square
%     drawnow
%     box on
%     ax = gca;
%     ax.LineWidth = 2;
%     ax.FontSize = 18;
 
end

% %% Plot S - Model and Experiment
% load('ExperimentAndModel_12PulseEnergies.mat');
% 
% open('Fig3b.fig');
% 
% tWZeroInd = find(tW>=0,1,'first');
% 
% for PulseEnergyInd = 3:Length
%     
%     s = S_Model_HR_All(PulseEnergyInd,:);
%     [~,sMaxInd] = max(s);
%     s = circshift(s,tWZeroInd-sMaxInd);
%     
%     figure(100)
%     hold on
%     plot(s + (PulseEnergyInd-1)-2,tW,'linewidth',2,'color',cmap(PulseEnergyInd-2,:))
%     hold off
%     ylabel('\Deltat [ps]')
%     xlabel('Energy loss [eV]')
%     ylim([-0.5 2])
%     xlim([-0.5 10])
%     ax = gca;
%     ax.FontSize = 18;
%     ax.LineWidth = 2;
%     ax.YDir = 'reverse';
%     box on
%     if mod(PulseEnergyInd,4) == 0
%     annotation('textbox',[0.17+((PulseEnergyInd-1)-2)/13.6, 0.17, 0.1, 0.1],'String',[num2str(PulseListNum(PulseEnergyInd)),' nJ'],...
%                 'color',cmap(PulseEnergyInd-2,:),'LineStyle','none','FontSize',18)
%     else
%     annotation('textbox',[0.17+((PulseEnergyInd-1)-2)/13.6, 0.17, 0.1, 0.1],'String',[num2str(PulseListNum(PulseEnergyInd)),],...
%                 'color',cmap(PulseEnergyInd-2,:),'LineStyle','none','FontSize',18)
%     end
% end
% 
%% Plot PV - Model and Experiment
load('Model_28PulseEnergies.mat');
open('Fig3c.fig');

Length = length(PulseEnergy);

for PulseEnergyInd = 1:Length
    PV_model(PulseEnergyInd) = max(S_Model_HR_All(PulseEnergyInd,:)) - min(S_Model_HR_All(PulseEnergyInd,:));
end

figure(1)
hold on
plot(PulseEnergy.*1e9./PulseEnergyGainFactor,PV_model,'-r','linewidth',2)
hold off
xlabel('Laser pulse energy [nJ]')
ylabel('P-V energy shift [eV]')
xlim([0 10.3])
ylim([0.2 1])
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 2;
box on
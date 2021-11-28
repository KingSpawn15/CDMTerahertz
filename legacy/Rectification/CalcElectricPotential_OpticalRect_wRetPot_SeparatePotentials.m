function [VCell] = CalcElectricPotential_OpticalRect_wRetPot_SeparatePotentials(PulseEnergy,t,z,x0,y0,Pol,v)

%% Physical constants
c = 299792458;%[m/s];
epsilon0 = 8.8541878128e-12;%[F/m]
mu0 = 1.25663706212e-6;%[H/m];
heta0 = 377;%[Ohm]

%% InAs parameters
d14 = 237.6e-12;%[V/m]
% Linear extrapolation in frequency, based on data from:
% Weber, Marvin J. Handbook of optical materials. CRC press, 2018.
% Assuming omega = 1.885e13[s^-1] or lambda = 100[um]

alpha = 7e6;%[m^-1]

%% Laser parameters
LaserSpotFWHM = 40e-6;%[m]
LaserSpotSigma = LaserSpotFWHM./(2*sqrt(2*log(2)));%[m]

LaerPulseTimeFWHM = 50e-15;%[s]
LaserPulseTimeSigma = LaerPulseTimeFWHM/(2*sqrt(2*log(2)));%[s]

%% Green's function variables (x',y',z') + kernel
xprime = (-1:5e-2:1)*round(LaserSpotSigma*3,5);%[m]
yprime = (0:5e-2:1)*1e-6;%[m]
zprime = xprime;%[m]

z0 = 30e-6;%[m] Sample +/-z boundary

[Xprime,Yprime,Zprime,Z] = ndgrid(xprime,yprime,zprime,z);

GreenKernel = ((x0-Xprime).^2+(y0-Yprime).^2+(Z-Zprime).^2).^(-1/2);

%% Down sample t to improve run speed
tC = t(t>-0.5e-12);
tCSubSampled = tC(1:30:end);

%% Calculate electric potential
y0Ind = find(yprime>=0,1,'first');
mz0Ind = find(zprime>=-z0,1,'first');
pz0Ind = find(zprime>=z0,1,'first');

LaserXZ = exp(-(Xprime.^2+Zprime.^2)./(LaserSpotSigma.^2));
LaserXZ = LaserXZ .* (Zprime<=z0) .* (Zprime>=-z0);

E0Squared = (PulseEnergy/LaerPulseTimeFWHM)*2*heta0;

t0 = 0.5e-12;%[s]
dt = t(2) - t(1);

VCell = cell(5,1);%zeros(length(tCSubSampled),length(z));

%parfor IndTime = 1:length(tCSubSampled)
for IndTime = 1:length(tCSubSampled)
    
    Tprime = tCSubSampled(IndTime) - (1/c)*sqrt((x0-Xprime).^2+(y0-Yprime).^2+(Z-Zprime).^2);

    LaserT = exp(-(Tprime-t0).^2./LaserPulseTimeSigma.^2);
    LaserT(Tprime<0) = 0;
    
    rho_xz = (1/sqrt(3)).*d14*E0Squared.*LaserT.*LaserXZ.*exp(-alpha.*Yprime).*...               
          ((2*sqrt(2)/LaserSpotSigma^2).*(Xprime.*sin(2*Pol)+Zprime.*cos(2*Pol)));
    rho_y = (1/sqrt(3)).*d14*E0Squared.*LaserT.*LaserXZ.*exp(-alpha.*Yprime).*...               
          ( - alpha);
      
    dPhi_xz = (1/(4*pi*epsilon0)).*rho_xz.*GreenKernel;
    dPhi_y = (1/(4*pi*epsilon0)).*rho_y.*GreenKernel;
    
    rho_Y0 = (1/sqrt(3)).*d14*E0Squared.*LaserT.*LaserXZ;
    rho_mZ0 = -d14*E0Squared.*LaserT.*LaserXZ.*exp(-alpha.*Yprime).*((2/sqrt(6)).*cos(2*Pol));      
    rho_pZ0 = -rho_mZ0;      
    
    dPhi_Y0 = (1/(4*pi*epsilon0)).*rho_Y0(:,y0Ind,:,:).*GreenKernel(:,y0Ind,:,:);
    dPhi_Z0 = (1/(4*pi*epsilon0)).*(rho_mZ0(:,:,mz0Ind,:).*GreenKernel(:,:,mz0Ind,:) + rho_pZ0(:,:,pz0Ind,:).*GreenKernel(:,:,pz0Ind,:));
    
    dPzdt = d14*E0Squared.*LaserT.*LaserXZ.*exp(-alpha.*Yprime).*...                
            (2/sqrt(6)).*cos(2*Pol).*(-2.*(Tprime-t0)./LaserPulseTimeSigma.^2);
    dPzdt(Tprime<0) = 0;
    
    dA = (mu0/(4*pi)).*dPzdt.*GreenKernel;

%     V(IndTime,:) = trapz(zprime,trapz(xprime, (trapz(yprime,(dPhi-v.*dA),2) + dPhi_Y0) ,1),3)...
%                    + trapz(xprime,trapz(yprime,dPhi_Z0,2),1);

    VCell{1}(IndTime,:) = trapz(zprime,trapz(yprime,trapz(xprime,dPhi_xz,1),2),3);
    VCell{2}(IndTime,:) = trapz(zprime,trapz(yprime,trapz(xprime,dPhi_y,1),2),3);
    VCell{3}(IndTime,:) = trapz(zprime,trapz(xprime,dPhi_Y0,1),3);
    VCell{4}(IndTime,:) = trapz(yprime,trapz(xprime,dPhi_Z0,1),2);
    VCell{5}(IndTime,:) = trapz(zprime,trapz(yprime,trapz(xprime,-v*dA,1),2),3);
 
end

%% Return to original t
for i = 1:5
    V = VCell{i};
    V = interp1(tCSubSampled,V,tC);
    V = V';
    V = [zeros(length(z),length(t)-length(tC)),V];
    V(isnan(V)) = 0;
    V = circshift(V,-round(t0/dt),2);
    VCell{i} = V;
end

end
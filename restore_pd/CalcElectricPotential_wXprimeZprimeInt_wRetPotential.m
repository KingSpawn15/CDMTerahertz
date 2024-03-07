function [Phi] = CalcElectricPotential_wXprimeZprimeInt_wRetPotential(PulseEnergy,t,z,x0,y0)

%% Physical constants
c = 299792458;%[m/s];
me0 = 9.109e-31;%[kg]
e = 1.60217662e-19;%[C]
epsilon0 = 8.8541878128e-12;%[F/m]

%% InAs parameters
alpha = 7e6;%[m^-1]
kappa = 12.3*epsilon0;
gamma = 3.3e12;%[s^-1]
me = 0.022*me0;%[kg]
mh = 0.6*me0;%[kg]
Lambda = 0.8;%[um]
hnew = 1.24./Lambda;%[eV]
Eg = 0.354;%[eV]
alphaGamma = 2.2;%[eV^-1]

Ee = (2.*(hnew-Eg)*mh)./(me+mh+sqrt((me+mh)^2+4*alphaGamma.*(hnew-Eg).*me.*mh));%[eV]
vte = sqrt((2.*Ee.*e.*(1+alphaGamma.*Ee))./(3.*me.*(1+4.*alphaGamma.*Ee.*(1+alphaGamma.*Ee))));%[m/s]
vt = vte;%[m/s]
me0tilda = me*3*(1+4*alphaGamma*Ee*(1+alphaGamma*Ee))^(3/2)/(3+8*alphaGamma*Ee*(1+alphaGamma*Ee));%[kg]

%% Doping
n_eq = 1e17*1e6;%[m^-3]%From resistivity measurement!
m_eq = mh;%p-type InAs

%% Laser parameters
LaserSpotFWHM = 40e-6;%[m]
LaserSpotSigma = LaserSpotFWHM./(2*sqrt(2*log(2)));%[m]
Alaser = (pi/4)*LaserSpotFWHM^2;%[m^2]
ExcitedVolume = Alaser*alpha^(-1);%[m^3]
n_exc = PulseEnergy/e/hnew/ExcitedVolume;%[m^-3]

%% Green's function variables (x',y',z') + kernel
xprime = (-1:4e-2:1)*LaserSpotFWHM*3/(2*sqrt(2*log(2)));%[m]
yprime = (0:4e-2:1)*1e-6;%[m]
zprime = xprime;%[m]

[Xprime,Yprime,Zprime,Z] = ndgrid(xprime,yprime,zprime,z);

GreenKernel = ((x0-Xprime).^2+(y0-Yprime).^2+(Z-Zprime).^2).^(-3/2);

%% Down sample t to improve run speed
tC = t(t>-0.5e-12);
tCSubSampled = tC(1:30:end);

%% Calculate electric potential
GaussianLaserSpot = exp(-(Xprime.^2+Zprime.^2)./(2*LaserSpotSigma.^2));

omega_y = sqrt( (e^2/kappa).*(n_exc.*exp(-alpha.*Yprime).*(1/me0tilda+1/mh)+n_eq/m_eq) - (gamma/2)^2 );

Phi = zeros(length(tCSubSampled),length(z));

parfor TimeInd = 1:length(tCSubSampled)
    
    Tprime = tCSubSampled(TimeInd) - (1/c)*sqrt((x0-Xprime).^2+(y0-Yprime).^2+(Z-Zprime).^2);
    
    g = (e*alpha*vt^2*n_exc.*exp(-alpha.*Yprime)./(omega_y.*(gamma^2+4.*omega_y.^2))).*...
                        (4.*omega_y.*cos(omega_y.*Tprime)+2*gamma.*sin(omega_y.*Tprime)).*...
                        exp(-gamma.*Tprime./2);
    g(Tprime<0) = 0;

    dPhi = GaussianLaserSpot.*GreenKernel.*g.*(y0-Yprime);  
    Phi(TimeInd,:) = (1/(4*pi*epsilon0)).* trapz(zprime,trapz(yprime,trapz(xprime,dPhi,1),2),3);           
    
end

%% Return to original t
Phi = interp1(tCSubSampled,Phi,tC);
Phi = Phi';
Phi = [zeros(length(z),length(t)-length(tC)),Phi];
Phi(isnan(Phi)) = 0;

end
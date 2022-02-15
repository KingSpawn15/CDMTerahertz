classdef IndiumArsenide < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        alpha
        kappa
        gamma
        me
        mh
        lambda
        hnew
        eg
        alpha_gamma
        n_eq
        m_eq
        epsilon_e
        v_t
        me0tilda
        d14
        gamma_factor
        phase

    end
    
    methods
        function self = IndiumArsenide()
            self.alpha = 7e6;%[m^-1]
            
            const = utils.constants_fundamantal();
            EPSILON_0 = const.('EPSILON_0');
            M_E = const.('M_E');
            self.kappa = 12.3*EPSILON_0;
            self.gamma = 1.3e13;%[s^-1]
            self.gamma_factor = 0.27;
            self.phase = -0.81;
            self.me = 0.022*M_E;%[kg]
            self.mh = 0.6*M_E;%[kg]
            self.lambda = 0.8;%[um]
            self.hnew = 1.24./self.lambda;%[eV]
            self.eg = 0.354;%[eV]
            self.alpha_gamma = 2.2;%[eV^-1]
            self.n_eq = 1e17*1e6;%[m^-3]%From resistivity measurement!
            self.m_eq = self.mh;%p-type InAs
            
            self.epsilon_e = self.photoexcited_electron_energy();
            self.v_t = self.velocity_t();
            self.me0tilda = self.photoelectron_mass();
            
            self.d14 = 237.6e-12;%[V/m]
            % Linear extrapolation in frequency, based on data from:
            % Weber, Marvin J. Handbook of optical materials. CRC press, 2018.
            % Assuming omega = 1.885e13[s^-1] or lambda = 100[um]
        end
        
        function v_t = velocity_t(self)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            const = utils.constants_fundamantal();
            Q_E = const.('Q_E');
            
            v_t = sqrt((2.*self.epsilon_e.*Q_E.*(1+self.alpha_gamma.*self.epsilon_e))...
                ./(3.*self.me.*(1+4.*self.alpha_gamma.*self.epsilon_e.*(1+self.alpha_gamma.*self.epsilon_e))));%[m/s]
            
           
        end
        
        
        function epsilon_e = photoexcited_electron_energy(self)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            epsilon_e = (2.*(self.hnew-self.eg)* self.mh)./...
                (self.me+self.mh+sqrt((self.me+self.mh)^2+4*self.alpha_gamma.*(self.hnew-self.eg).*self.me.*self.mh));%[eV]
        end


        function me0tilda = photoelectron_mass(self)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
                                 
            me0tilda = self.me*3*(1+4*self.alpha_gamma*self.epsilon_e*(1+self.alpha_gamma*self.epsilon_e))^(3/2)...
                /(3+8*self.alpha_gamma*self.epsilon_e*(1+self.alpha_gamma*self.epsilon_e));%[kg]
            
        end
        
        function omega_y = calculate_omega_y(self, n_exc, yprime_grid, gaussian_laser_spot)
            const = utils.constants_fundamantal();
            Q_E = const.('Q_E');
            
            omega_y = sqrt( (Q_E^2/self.kappa).*(n_exc.*gaussian_laser_spot.*exp(-self.alpha.*yprime_grid).*...
                (1/self.me0tilda+1/self.mh)+self.n_eq/self.m_eq) - (self.gamma/2)^2 );
        end
        
    end
end


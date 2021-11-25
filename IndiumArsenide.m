classdef IndiumArsenide
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

    end
    
    methods
        function self = IndiumArsenide()
            self.alpha = 7e6;%[m^-1]
            self.kappa = 12.3*epsilon0;
            self.gamma = 3.3e12;%[s^-1]
            self.me = 0.022*me0;%[kg]
            self.mh = 0.6*me0;%[kg]
            self.lambda = 0.8;%[um]
            self.hnew = 1.24./lambda;%[eV]
            self.eg = 0.354;%[eV]
            self.alpha_gamma = 2.2;%[eV^-1]
            self.n_eq = 1e17*1e6;%[m^-3]%From resistivity measurement!
            self.m_eq = mh;%p-type InAs
            
            self.epsilon_e = self.photoexcited_electron_energy(self);
            self.v_t = self.velocity_t();
            self.me0tilda = photoelectron_mass();
        end
        
        function v_t = velocity_t(self)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            
            v_t = sqrt((2.*self.epsilon_e.*e.*(1+self.alpha_gamma.*self.epsilon_e))...
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
        
    end
end


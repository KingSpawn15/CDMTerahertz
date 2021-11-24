classdef DiscretizationSteps
    %DISCRETIZATIONSTEPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ddt
        deltat
        z
        deltaz
        t
        omega
        energy
    end
    
    methods
        
        function self = DiscretizationSteps(ddt, delay_max, ...
                ddz, zmax,...
                fs, l)
            [self.ddt, self.deltat] = self.delay(ddt, delay_max);
            [self.z , self.deltaz] = self.zstep(ddz, zmax);
            [self.t, self.omega, self.energy] = self.incidence(fs, l);
        end

    end

    methods(Static)
        function [ddt, deltat] = delay(ddt, delay_max)
            %DISCRETIZATION_DT Summary of this function goes here
            %   Detailed explanation goes here
            
            
            % Time delay
            
            deltat = -delay_max:ddt:delay_max;
            
            
        end
        
        function [z , deltaz] = zstep(ddz, zmax)
            %DISCRETIZATION_Z Summary of this function goes here
            %   Detailed explanation goes here
            
            
            z = -zmax : ddz : zmax;%[m]
            deltaz = z(2)-z(1);
            
            
        end
        
        function [t, omega, energy] = incidence(fs, l)
            
            const = utils.constants_fundamantal();
            HBAR = const.('HBAR');
            Q_E = const.('Q_E');
            dt = 1/fs;%[s]
            domega = fs/l*(2*pi);%[Hz]
            t = (-l/2:l/2).*dt;%[s]
            omega = (-l/2:l/2)*domega;%[Hz]
            energy = omega.*HBAR/Q_E;%[eV]
            
        end
        
    end
    
end
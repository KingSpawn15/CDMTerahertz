classdef Discretization < handle
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
        x0
        y0
        t0
        xprime
        yprime
        zprime
        d_xprime
        d_yprime
        d_zprime
        XPRIME
        YPRIME
        ZPRIME
        Z
        z_max
        dt
    end
    
    methods
        
        function self = Discretization(discretization_params)
            
            self.x0 = discretization_params.x0;
            self.y0 = discretization_params.y0;
            self.d_xprime = discretization_params.d_xprime;
            self.d_yprime = discretization_params.d_yprime;
            self.d_zprime = discretization_params.d_zprime;
            [self.ddt, self.deltat] = self.delay(discretization_params.ddt, discretization_params.delay_max);
            [self.z , self.deltaz] = self.zstep(discretization_params.ddz, discretization_params.zmax);
            [self.t, self.omega, self.energy] = self.incidence(discretization_params.fs, discretization_params.l); 
            self.dt = self.t(2) - self.t(1);
            self.t0 = discretization_params.t0;
            [self.xprime , self.yprime , self.zprime , ...
                self.XPRIME, self.YPRIME, self.ZPRIME, self.Z] = ...
                self.space_points_grid(discretization_params , self.z);
            self.z_max = discretization_params.z_max;
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
    
    methods(Static)
        function [xprime , yprime , zprime , ...
                XPRIME, YPRIME, ZPRIME, Z] = space_points_grid(discretization_params , z)
            
            xprime = - discretization_params.xprime_max : discretization_params.d_xprime :  discretization_params.xprime_max ;
            yprime = 0 : discretization_params.d_yprime :  discretization_params.yprime_max ;
            zprime = - discretization_params.zprime_max : discretization_params.d_zprime :  discretization_params.zprime_max ;
            
            [XPRIME,YPRIME,...
                ZPRIME,Z] = ...
                ndgrid(xprime,yprime,zprime,...
                z);
            
        end
    end
end
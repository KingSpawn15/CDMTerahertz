classdef Laser
    %LASER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pulse_energy_gain_factor
        theta_pol
        pulse_energy
        laser_spot_fwhm
        laser_spot_sigma
    end
    
    methods
        
        function self = Laser(laser_parameters)
            
            self.pulse_energy = laser_parameters.pulse_energy_gain_factor * laser_parameters.pulse_energy_experiment;
            self.theta_pol = laser_parameters.theta_pol;
            
            if isfield(laser_parameters , 'laser_spot_fwhm')
                self.laser_spot_fwhm = laser_parameters.laser_spot_fwhm;
                self.laser_spot_sigma = self.calculate_sigma(self.laser_spot_fwhm);
            end
        end
        
    end
    
    methods(Static)
        
        function sigma = calculate_sigma(fwhm)
            sigma = fwhm./(2*sqrt(2*log(2)));
        end
        
    end
    
end


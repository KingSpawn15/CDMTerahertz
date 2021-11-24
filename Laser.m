classdef Laser
    %LASER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pulse_energy_gain_factor
        theta_pol
        pulse_energy
        
    end
    
    methods
        
        function self = Laser(pulse_energy_experiment, pulse_energy_gain_factor, theta_pol)
            
            self.pulse_energy = pulse_energy_gain_factor * pulse_energy_experiment;
            self.theta_pol = theta_pol;
            
        end
        
    end
end


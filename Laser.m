classdef Laser
    %LASER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pulse_energy_gain_factor
        theta_pol
        pulse_energy
        laser_spot_fwhm
        laser_spot_sigma
        laser_pulse_time_fwhm
        laser_pulse_time_sigma
    end
    
    methods
        
        function self = Laser(laser_parameters)
            
            self.pulse_energy = laser_parameters.pulse_energy_gain_factor * laser_parameters.pulse_energy_experiment;
            self.theta_pol = laser_parameters.theta_pol;
            
            if isfield(laser_parameters , 'laser_spot_fwhm')
                self.laser_spot_fwhm = laser_parameters.laser_spot_fwhm;
                self.laser_spot_sigma = self.calculate_sigma(self.laser_spot_fwhm);
            end
            
            if isfield(laser_parameters , 'laser_pulse_time_fwhm')
                self.laser_pulse_time_fwhm = laser_parameters.laser_pulse_time_fwhm;%[s]
                self.laser_pulse_time_sigma = self.calculate_sigma(laser_parameters.laser_pulse_time_fwhm);%[s]
            end
            
        end
       
        function n_exc = excited_carriers(self, alpha, hnew)
            
            const = utils.constants_fundamantal();
            Q_E = const.('Q_E');
            excited_volume = (pi/4)*self.laser_spot_fwhm^2*alpha^(-1);%[m^3]
            n_exc = self.pulse_energy/Q_E/hnew/excited_volume;%[m^-3]

        end

    end
    
    methods(Static)
        
        function sigma = calculate_sigma(fwhm)
            sigma = fwhm./(2*sqrt(2*log(2)));
        end
        
    end
    
end


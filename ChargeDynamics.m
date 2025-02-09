classdef ChargeDynamics
    %CDEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        
        
        function interaction_v =  ...
                interaction_potential_rectification(discretization, material,...
                laser , electron, subsampling)    
            % Physical constants
            const = utils.constants_fundamantal();
            C = const.('C');
            EPSILON_0 = const.('EPSILON_0');
            MU_0 = const.('MU_0');
            ETA_0 = const.('ETA_0');
            
            % InAs parameters
            d14 = material.d14;
            alpha = material.alpha;
            
            % Laser parameters
            
            laser_spot_sigma = laser.laser_spot_sigma;
            
            laser_pulse_time_fwhm = laser.laser_pulse_time_fwhm;%[s]
            laser_pulse_time_sigma = laser.laser_pulse_time_sigma;%[s]
            pulse_energy = laser.pulse_energy;
            
            % Green's function variables (x',y',z') + kernel
            %             xprime = discretization.xprime;
            yprime = discretization.yprime;%[m]
            zprime = discretization.zprime;%[m]
            
            d_xprime = discretization.d_xprime;
            d_yprime = discretization.d_yprime;
            d_zprime = discretization.d_xprime;
            
            XPRIME = discretization.XPRIME;
            YPRIME = discretization.YPRIME;
            ZPRIME= discretization.ZPRIME;
            
            Z = discretization.Z;
            z_max = discretization.z_max;%[m] Sample +/-z boundary
            
            green_kernel  = ChargeDynamics.calculate_green_kernel_rectification(discretization);
            
            % Down sample t to improve run speed
            
            t_c = discretization.t(discretization.t > discretization.t0);
            t_c_subsampled = t_c(1:subsampling.tc_subsampling:end);
            
            
            % Calculate electric potential
            y0Ind = find(yprime>=0,1,'first');
            mz0Ind = find(zprime>=-z_max,1,'first');
            pz0Ind = find(zprime>=z_max,1,'first');
            
            LaserXZ = exp(-(XPRIME.^2+ZPRIME.^2)./(laser_spot_sigma.^2));
            LaserXZ = LaserXZ .* (ZPRIME<=z_max) .* (ZPRIME>=-z_max);
            
            E0Squared = (pulse_energy/laser_pulse_time_fwhm)*2*ETA_0;
            
            t0 = -discretization.t0;%[s]  %%% WTF is this line
            dt = discretization.dt;
            
            interaction_v = zeros(length(t_c_subsampled),length(discretization.z));
            
            theta_pol = laser.theta_pol;
            electron_velocity = electron.electron_velocity;
            t_r = (1/C)*sqrt((discretization.x0-XPRIME).^2+(discretization.y0-YPRIME).^2+(Z-ZPRIME).^2);
            parfor time_ind = 1:length(t_c_subsampled)
                
                disp(strcat('time_ind', num2str(time_ind), ' out of ', num2str(length(t_c_subsampled))));
                t_prime = t_c_subsampled(time_ind) - t_r;
                
                laser_t = exp(-(t_prime-t0).^2./laser_pulse_time_sigma.^2);
                laser_t(t_prime<0) = 0;
                
                rho = (1/sqrt(3)).*d14*E0Squared.*laser_t.*LaserXZ.*exp(-alpha.*YPRIME).*...
                    ((2*sqrt(2)/laser_spot_sigma^2).*(XPRIME.*sin(2*theta_pol)+ZPRIME.*cos(2*theta_pol)) - alpha);
                
                dPhi = (1/(4*pi*EPSILON_0)).*rho.*green_kernel;
                
                rho_Y0 = (1/sqrt(3)).*d14*E0Squared.*laser_t.*LaserXZ;
                rho_mZ0 = -d14*E0Squared.*laser_t.*LaserXZ.*exp(-alpha.*YPRIME).*((2/sqrt(6)).*cos(2*theta_pol));
                rho_pZ0 = -rho_mZ0;
                
                dPhi_Y0 = (1/(4*pi*EPSILON_0)).*rho_Y0(:,y0Ind,:,:).*green_kernel(:,y0Ind,:,:);
                dPhi_Z0 = (1/(4*pi*EPSILON_0)).*(rho_mZ0(:,:,mz0Ind,:).*green_kernel(:,:,mz0Ind,:) + rho_pZ0(:,:,pz0Ind,:).*green_kernel(:,:,pz0Ind,:));
                
                dPzdt = d14*E0Squared.*laser_t.*LaserXZ.*exp(-alpha.*YPRIME).*...
                    (2/sqrt(6)).*cos(2*theta_pol).*(-2.*(t_prime-t0)./laser_pulse_time_sigma.^2);
                dPzdt(t_prime<0) = 0;
                
                dA = (MU_0/(4*pi)).*dPzdt.*green_kernel;
                
                %     V(time_ind,:) = trapz(zprime,trapz(xprime, (trapz(yprime,(dPhi-v.*dA),2) + dPhi_Y0) ,1),3)...
                %         + trapz(xprime,trapz(yprime,dPhi_Z0,2),1);
                if (time_ind == 207)
                    disp("");
                end
                interaction_v(time_ind,:) = trapz(d_zprime,trapz(d_xprime, (trapz(d_yprime,(dPhi-electron_velocity.*dA),2) + dPhi_Y0) ,1),3)...
                    + trapz(d_xprime,trapz(d_yprime,dPhi_Z0,2),1);
            end
            
            % Return to original t
            interaction_v = interp1(t_c_subsampled,interaction_v,t_c);
            interaction_v = interaction_v';
            interaction_v = [zeros(length(discretization.z),length(discretization.t)-length(t_c)),interaction_v];
            interaction_v(isnan(interaction_v)) = 0;
            interaction_v = circshift(interaction_v,-round(t0/dt),2);
            
        end
        
        
        function interaction_v =  ...
                interaction_potential_photodember(discretization, material,...
                laser , subsampling)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            const = utils.constants_fundamantal();
            EPSILON_0 = const.('EPSILON_0');
            C = const.('C');
            Q_E = const.('Q_E');
            
            t_c = discretization.t(discretization.t > discretization.t0);
            t_c_subsampled = t_c(1:subsampling.tc_subsampling:end);
            
            x0 = discretization.x0;
            y0 = discretization.y0;
            z = discretization.z;
            
            xprime = discretization.xprime;
            yprime = discretization.yprime;
            zprime = discretization.zprime;
            
            XPRIME = discretization.XPRIME;
            YPRIME = discretization.YPRIME;
            ZPRIME= discretization.ZPRIME;
            
            Z = discretization.Z;
            
            interaction_v =  zeros(length(t_c_subsampled),length(discretization.z));
            
            green_kernel  = ChargeDynamics.calculate_green_kernel(discretization);
            gaussian_laser_spot = exp(-(XPRIME.^2+ZPRIME.^2)./(2*laser.laser_spot_sigma.^2));
            
            l_tc = length(t_c_subsampled);
            
            n_exc = laser.excited_carriers(material.alpha, material.hnew);
            omega_y = material.calculate_omega_y(n_exc, YPRIME);
            
            t_r = (1/C)*sqrt((x0-XPRIME).^2+(y0-YPRIME).^2+(Z-ZPRIME).^2);
            
            alpha = material.alpha;
            gamma = material.gamma;
            v_t = material.v_t;
            
            parfor time_ind = 1:l_tc
                
                disp(strcat('time_ind', num2str(time_ind), ' out of ', num2str(l_tc)));
                t_prime = t_c_subsampled(time_ind) - t_r;
                
                g = (Q_E*alpha*v_t^2*...
                    n_exc.*exp(-alpha.*YPRIME)./(omega_y.*(gamma^2+4.*omega_y.^2))).*...
                    (4.*omega_y.*cos(omega_y.*t_prime)+2*gamma.*sin(omega_y.*t_prime)).*...
                    exp(-gamma.*t_prime./2);
                g(t_prime<0) = 0;
                
                dv = gaussian_laser_spot.*green_kernel.*g.*(y0-YPRIME);
                
                interaction_v(time_ind,:) = (1/(4*pi*EPSILON_0)).* trapz(zprime,trapz(yprime,trapz(xprime,dv,1),2),3);
                
            end
            
            
            interaction_v =  interp1(t_c_subsampled,interaction_v,t_c);
            interaction_v =  interaction_v';
            interaction_v =  [zeros(length(z),length(discretization.t)-length(t_c)),interaction_v];
            interaction_v(isnan(interaction_v)) = 0;
            
        end
        
        function green_kernel = calculate_green_kernel(discretization)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            green_kernel = ((discretization.x0-discretization.XPRIME).^2 + ...
                (discretization.y0-discretization.YPRIME).^2+(discretization.Z-discretization.ZPRIME).^2).^(-3/2);
            
        end
        
        function green_kernel = calculate_green_kernel_rectification(discretization)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            green_kernel = ((discretization.x0-discretization.XPRIME).^2 + ...
                (discretization.y0-discretization.YPRIME).^2+(discretization.Z-discretization.ZPRIME).^2).^(-1/2);
            
        end
        
        function g = dpdt(laser , material, t_prime, YPRIME)
            
            const = utils.constants_fundamantal();
            Q_E = const.('Q_E');
            
            n_exc = laser.excited_carriers(material.alpha, material.hnew);
            omega_y = material.calculate_omega_y(n_exc, YPRIME);
            
            
            g = (Q_E*material.alpha*material.v_t^2*...
                n_exc.*exp(-material.alpha.*YPRIME)./(omega_y.*(material.gamma^2+4.*omega_y.^2))).*...
                (4.*omega_y.*cos(omega_y.*t_prime)+2*material.gamma.*sin(omega_y.*t_prime)).*...
                exp(-material.gamma.*t_prime./2);
            g(t_prime<0) = 0;
            
        end
    end
end


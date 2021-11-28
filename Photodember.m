% classdef Photodember
%     %PHOTODEMBER Summary of this class goes here
%     %   Detailed explanation goes here
%     
%     properties
%         
%     end
%     
%     methods
%         function self = Photodember()
%             %PHOTODEMBER Construct an instance of this class
%             %   Detailed explanation goes here
%             
%         end
%         
%     end
%     
%     methods(Static)
%         
%         
%         function v = interaction_potential(discretization, material, laser , subsampling)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             
%             const = utils.constants_fundamantal();
%             EPSILON_0 = const.('EPSILON_0');
%             C = const.('C');
%             
%             t_c = discretization.t(discretization.t > discretization.t0);
%             t_c_subsampled = t_c(1:subsampling.tc_subsampling:end);
%             
%             v = zeros(length(t_c_subsampled),length(discretization.z));
%             
%              
%             green_kernel = Photodember.calculate_green_kernel(discretization);
%             gaussian_laser_spot = exp(-(discretization.XPRIME.^2+discretization.ZPRIME.^2)./(2*laser.laser_spot_sigma.^2));
%             
%             
%             l_tc = length(t_c_subsampled);
%             
%             x0 = discretization.x0;
%             y0 = discretization.y0;
%             Z = discretization.ZPRIME;
%             
%             zprime = discretization.zprime;
%             yprime = discretization.yprime;
%             xprime = discretization.xprime;
%             
%             XPRIME = discretization.XPRIME;
%             YPRIME = discretization.YPRIME;
%             ZPRIME = discretization.ZPRIME;
%             
%             const = utils.constants_fundamantal();
%             Q_E = const.('Q_E');
%             n_exc = laser.excited_carriers(material.alpha, material.hnew);
%             omega_y = material.calculate_omega_y(n_exc, YPRIME);
%             
%             t_r = (1/C)*sqrt((x0-XPRIME).^2+(y0-YPRIME).^2+(Z-ZPRIME).^2);
%             
%             parfor time_ind = 1:l_tc
%                 
%                 disp(strcat('IndTime', num2str(time_ind), ' out of ', num2str(l_tc)));
%                 t_prime = t_c_subsampled(time_ind) - t_r;
%                 %                 g = Photodember.dpdt(laser , material, t_prime, YPRIME);
%                 
%                 
%                 
%                 g = (Q_E*material.alpha*material.v_t^2*...
%                     n_exc.*exp(-material.alpha.*YPRIME)./(omega_y.*(material.gamma^2+4.*omega_y.^2))).*...
%                     (4.*omega_y.*cos(omega_y.*t_prime)+2*material.gamma.*sin(omega_y.*t_prime)).*...
%                     exp(-material.gamma.*t_prime./2);
%                 g(t_prime<0) = 0;
%                 
%                 
%                 dv = gaussian_laser_spot.*green_kernel.*g.*(y0-YPRIME);
%                 
%                 v(time_ind,:) = (1/(4*pi*EPSILON_0)).* trapz(zprime,trapz(yprime,trapz(xprime,dv,1),2),3);
%                 
%             end
%             
%             %% Return to original t
%             v = interp1(self.tc_subsampled,v,t_c);
%             v = v';
%             v = [zeros(length(z),length(t)-length(t_c)),v];
%             v(isnan(v)) = 0;
%             
%         end
%         
%         
%         function interaction_v =  interaction_potential_2(discretization, material, laser , subsampling)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             
%             const = utils.constants_fundamantal();
%             EPSILON_0 = const.('EPSILON_0');
%             C = const.('C');
%             Q_E = const.('Q_E');
%             
%             t_c = discretization.t(discretization.t > discretization.t0);
%             t_c_subsampled = t_c(1:subsampling.tc_subsampling:end);
%             
%             x0 = discretization.x0;
%             y0 = discretization.y0;
%             z = discretization.z;
%             
%             xprime = discretization.xprime;
%             yprime = discretization.yprime;
%             zprime = discretization.zprime;
%             
%             XPRIME = discretization.XPRIME;
%             YPRIME = discretization.YPRIME;
%             ZPRIME= discretization.ZPRIME;
%             
%             Z = discretization.Z;
%             
%             interaction_v =  zeros(length(t_c_subsampled),length(discretization.z));
%                          
%             green_kernel  = Photodember.calculate_green_kernel(discretization);
%             gaussian_laser_spot = exp(-(XPRIME.^2+ZPRIME.^2)./(2*laser.laser_spot_sigma.^2));
% 
%             l_tc = length(t_c_subsampled);
%             
%             n_exc = laser.excited_carriers(material.alpha, material.hnew);
%             omega_y = material.calculate_omega_y(n_exc, YPRIME);
%             
%             t_r = (1/C)*sqrt((x0-XPRIME).^2+(y0-YPRIME).^2+(Z-ZPRIME).^2);
%             
%             alpha = material.alpha;
%             gamma = material.gamma;
%             v_t = material.v_t;
%             
%             parfor time_ind = 1:l_tc
%                 
%                 disp(strcat('IndTime', num2str(time_ind), ' out of ', num2str(l_tc)));
%                 t_prime = t_c_subsampled(time_ind) - t_r;                
%                 
%                 g = (Q_E*alpha*v_t^2*...
%                     n_exc.*exp(-alpha.*YPRIME)./(omega_y.*(gamma^2+4.*omega_y.^2))).*...
%                     (4.*omega_y.*cos(omega_y.*t_prime)+2*gamma.*sin(omega_y.*t_prime)).*...
%                     exp(-gamma.*t_prime./2);
%                 g(t_prime<0) = 0;
%                 
%                 dv = gaussian_laser_spot.*green_kernel.*g.*(y0-YPRIME);
%                 
%                 interaction_v(time_ind,:) = (1/(4*pi*EPSILON_0)).* trapz(zprime,trapz(yprime,trapz(xprime,dv,1),2),3);
%                 
%             end
%             
%             %% Return to original t
%             interaction_v =  interp1(t_c_subsampled,interaction_v,t_c);
%             interaction_v =  interaction_v';
%             interaction_v =  [zeros(length(z),length(discretization.t)-length(t_c)),interaction_v];
%             interaction_v(isnan(interaction_v)) = 0;
%             
%         end
%         
%         function green_kernel = calculate_green_kernel(discretization)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             green_kernel = ((discretization.x0-discretization.XPRIME).^2 + ...
%                 (discretization.y0-discretization.YPRIME).^2+(discretization.Z-discretization.ZPRIME).^2).^(-3/2);
%             
%         end
%         
%         
%         function g = dpdt(laser , material, t_prime, YPRIME)
%             
%             const = utils.constants_fundamantal();
%             Q_E = const.('Q_E');
%             
%             n_exc = laser.excited_carriers(material.alpha, material.hnew);
%             omega_y = material.calculate_omega_y(n_exc, YPRIME);
%             
%             
%             g = (Q_E*material.alpha*material.v_t^2*...
%                 n_exc.*exp(-material.alpha.*YPRIME)./(omega_y.*(material.gamma^2+4.*omega_y.^2))).*...
%                 (4.*omega_y.*cos(omega_y.*t_prime)+2*material.gamma.*sin(omega_y.*t_prime)).*...
%                 exp(-material.gamma.*t_prime./2);
%             g(t_prime<0) = 0;
%             
%         end
%         
%     end
%     
%     
% end
% 

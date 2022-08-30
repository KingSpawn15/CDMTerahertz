classdef ChargeDynamics
    %CDEM Summary of this class goes here
    %   Detailed explanation goes here

    properties

    end

    methods(Static)



        function interaction_v =  ...
                interaction_potential_rectification(discretization, material,...
                laser , electron, numerical_parameters)
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
            t_c_subsampled = t_c(1:numerical_parameters.tc_subsampling:end);


            % Calculate electric potential
            y0_ind = find(yprime>=0,1,'first');
            mz0_ind = find(zprime>=-z_max,1,'first');
            pz0_ind = find(zprime>=z_max,1,'first');

            laser_xz = exp(-(XPRIME.^2+ZPRIME.^2)./(laser_spot_sigma.^2));
            laser_xz = laser_xz .* (ZPRIME<=z_max) .* (ZPRIME>=-z_max);

            e0_squared = (pulse_energy/laser_pulse_time_fwhm)*2*ETA_0;

            t0 = -discretization.t0;%[s]  %%% WTF is this line
            dt = discretization.dt;

            interaction_v = zeros(length(t_c_subsampled),length(discretization.z));

            theta_pol = laser.theta_pol;
            electron_velocity = electron.electron_velocity;
            t_r = (1/C)*sqrt((discretization.x0-XPRIME).^2+(discretization.y0-YPRIME).^2+(Z-ZPRIME).^2);

            disp("Calculating potential rectification ....");
            begin = '['; last = '] ... Rectification'; arrow = '>';
            body = '='; empty = ' ';
            parfor time_ind = 1:length(t_c_subsampled)

                body_disp = repmat(body , [1,fix(time_ind / length(t_c_subsampled) * 70)]);
                empty_disp = repmat(empty , [1,70 - fix(time_ind / length(t_c_subsampled) * 70)]);
                disp([begin, body_disp, arrow, empty_disp , last]);
                %                 disp([num2str(time_ind),'/',num2str(length(t_c_subsampled))]);

                t_prime = t_c_subsampled(time_ind) - t_r;

                laser_t = exp(-(t_prime-t0).^2./laser_pulse_time_sigma.^2);
                sigma_side = 2.5 * laser_pulse_time_sigma;
                amplitude_side = 0.005;
                t0_side = 8* 0.6 * sigma_side;
               
                laser_t = laser_t + ...
                    1* amplitude_side.*exp(-(t_prime - t0 -t0_side).^2./(sigma_side).^2) + ...
                     1 * amplitude_side.*exp(-(t_prime - t0 +t0_side).^2./(sigma_side).^2);
                %                 laser_t = laser_t .* 1./(1 + exp(+(t_prime-t0)));
                laser_t(t_prime<0) = 0;

                lp = -laser_pulse_time_sigma;
                lp = (1/laser_pulse_time_fwhm);

                %                 asymmetry_pulse = 1./(1 + exp(lp.*(t0 - ...
                %                     t_prime + laser_pulse_time_fwhm*0)));
                %                 laser_t_asym = laser_t .* asymmetry_pulse;

                b = 0;
                asymmetry_pulse = 1./(1 + b*(t_prime-t0)./sqrt((t_prime-t0).^2 + laser_pulse_time_fwhm.^2));
                laser_t_asym_2 = laser_t .^ asymmetry_pulse;

%                 laser_t_asym_2 = laser_t;
                laser_t_asym = laser_t_asym_2;


                rho = (1/sqrt(3)).*d14*e0_squared.*laser_t_asym.*laser_xz.*exp(-alpha.*YPRIME).*...
                    ((2*sqrt(2)/laser_spot_sigma^2).*(XPRIME.*sin(2*theta_pol)+ZPRIME.*cos(2*theta_pol)) - alpha);



                dPhi = .02*(1/(4*pi*EPSILON_0)).*rho.*green_kernel;

                rho_Y0 = 0*(1/sqrt(3)).*d14*e0_squared.*laser_t_asym.*laser_xz;
                rho_mZ0 = -d14*e0_squared.*laser_t_asym.*laser_xz.*exp(-alpha.*YPRIME).*((2/sqrt(6)).*cos(2*theta_pol));
                rho_pZ0 = -rho_mZ0;

                dPhi_Y0 = (1/(4*pi*EPSILON_0)).*rho_Y0(:,y0_ind,:,:).*green_kernel(:,y0_ind,:,:);
                dPhi_Z0 = 50*(1/(4*pi*EPSILON_0)).*(rho_mZ0(:,:,mz0_ind,:).*green_kernel(:,:,mz0_ind,:) + rho_pZ0(:,:,pz0_ind,:).*green_kernel(:,:,pz0_ind,:));
                %
%                 dPzdt = d14*e0_squared.*laser_t.*laser_xz.*exp(-alpha.*YPRIME).*...
%                     (2/sqrt(6)).*cos(2*theta_pol).*(-2.*(t_prime-t0)./laser_pulse_time_sigma.^2);

                dPzdt = d14*e0_squared.*laser_t_asym_2.*laser_xz.*exp(-alpha.*YPRIME).*...
                    (2/sqrt(6)).*cos(2*theta_pol).*(-2.*(t_prime-t0)./laser_pulse_time_sigma.^2);


                %                 dPzdt = d14*e0_squared.*laser_t.*laser_xz.*exp(-alpha.*YPRIME).*...
                %                     (2/sqrt(6)).*cos(2*theta_pol).*((exp(lp.*(t0 - t_prime)).*(laser_pulse_time_fwhm.^2.*lp + 2.*t0 - ...
                %                     2.*t_prime) + 2.*(t0 - t_prime))./((1 + exp(lp.*(t0 - ...
                %                     t_prime))).^2.*laser_pulse_time_fwhm.^2));

%                 dPzdt(t_prime<0) = 0;

                dA = 1 * (MU_0/(4*pi)).*dPzdt.*green_kernel;

                
%                                     interaction_v(time_ind,:) = trapz(d_zprime,trapz(d_xprime, (trapz(d_yprime,(dPhi-electron_velocity.*dA),2) + dPhi_Y0) ,1),3)...
%                                         + trapz(d_xprime,trapz(d_yprime,dPhi_Z0,2),1);

                interaction_v(time_ind,:) = trapz(d_zprime,trapz(d_xprime, (trapz(d_yprime,(-electron_velocity.*dA),2)) ,1),3);
            end

            % Return to original t
            interaction_v = interp1(t_c_subsampled,interaction_v,t_c);
            interaction_v = interaction_v';
            interaction_v = [zeros(length(discretization.z),length(discretization.t)-length(t_c)),interaction_v];
            interaction_v(isnan(interaction_v)) = 0;
            interaction_v = circshift(interaction_v,-round(t0/dt),2);

        end



        function interaction_v =  ...
                interaction_potential_rectification_mod(discretization, material,...
                laser , electron, numerical_parameters)
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
            t_c_subsampled = t_c(1:numerical_parameters.tc_subsampling:end);


            % Calculate electric potential
            y0_ind = find(yprime>=0,1,'first');
            mz0_ind = find(zprime>=-z_max,1,'first');
            pz0_ind = find(zprime>=z_max,1,'first');

            laser_xz = exp(-(XPRIME.^2+ZPRIME.^2)./(laser_spot_sigma.^2));
            laser_xz = laser_xz .* (ZPRIME<=z_max) .* (ZPRIME>=-z_max);

            e0_squared = (pulse_energy/laser_pulse_time_fwhm)*2*ETA_0;

            t0 = -discretization.t0;%[s]  %%% WTF is this line
            dt = discretization.dt;

            interaction_v = zeros(length(t_c_subsampled),length(discretization.z));

            theta_pol = laser.theta_pol;
            electron_velocity = electron.electron_velocity;
            t_r = (1/C)*sqrt((discretization.x0-XPRIME).^2+(discretization.y0-YPRIME).^2+(Z-ZPRIME).^2);

            disp("Calculating potential rectification ....");
            begin = '['; last = '] ... Rectification'; arrow = '>';
            body = '='; empty = ' ';
            parfor time_ind = 1:length(t_c_subsampled)

                body_disp = repmat(body , [1,fix(time_ind / length(t_c_subsampled) * 70)]);
                empty_disp = repmat(empty , [1,70 - fix(time_ind / length(t_c_subsampled) * 70)]);
                disp([begin, body_disp, arrow, empty_disp , last]);
                %                 disp([num2str(time_ind),'/',num2str(length(t_c_subsampled))]);

                t_prime = t_c_subsampled(time_ind) - t_r;

                laser_t = exp(-(t_prime-t0).^2./laser_pulse_time_sigma.^2);
                %                 laser_t = laser_t .* 1./(1 + exp(+(t_prime-t0)));
                laser_t(t_prime<0) = 0;



                rho = (1/sqrt(3)).*d14*e0_squared.*laser_t_asym.*laser_xz.*exp(-alpha.*YPRIME).*...
                    ((2*sqrt(2)/laser_spot_sigma^2).*(XPRIME.*sin(2*theta_pol)+ZPRIME.*cos(2*theta_pol)) - alpha);



                dPhi = (1/(4*pi*EPSILON_0)).*rho.*green_kernel;

                rho_Y0 = (1/sqrt(3)).*d14*e0_squared.*laser_t_asym.*laser_xz;
                rho_mZ0 = -d14*e0_squared.*laser_t_asym.*laser_xz.*exp(-alpha.*YPRIME).*((2/sqrt(6)).*cos(2*theta_pol));
                rho_pZ0 = -rho_mZ0;

                dPhi_Y0 = (1/(4*pi*EPSILON_0)).*rho_Y0(:,y0_ind,:,:).*green_kernel(:,y0_ind,:,:);
                dPhi_Z0 = (1/(4*pi*EPSILON_0)).*(rho_mZ0(:,:,mz0_ind,:).*green_kernel(:,:,mz0_ind,:) + rho_pZ0(:,:,pz0_ind,:).*green_kernel(:,:,pz0_ind,:));
                %
                %                                 dPzdt = d14*e0_squared.*laser_t.*laser_xz.*exp(-alpha.*YPRIME).*...
                %                                     (2/sqrt(6)).*cos(2*theta_pol).*(-2.*(t_prime-t0)./laser_pulse_time_sigma.^2);
                %
                %                 dPzdt = d14*e0_squared.*laser_t_asym_2.*laser_xz.*exp(-alpha.*YPRIME).*...
                %                     (2/sqrt(6)).*cos(2*theta_pol).*(-2.*(t_prime-t0)./laser_pulse_time_sigma.^2);
                %

                dPzdt = d14*e0_squared.*laser_t.*laser_xz.*exp(-alpha.*YPRIME).*...
                    (2/sqrt(6)).*cos(2*theta_pol).*((exp(lp.*(t0 - t_prime)).*(laser_pulse_time_fwhm.^2.*lp + 2.*t0 - ...
                    2.*t_prime) + 2.*(t0 - t_prime))./((1 + exp(lp.*(t0 - ...
                    t_prime))).^2.*laser_pulse_time_fwhm.^2));

                dPzdt(t_prime<0) = 0;

                dA = -1 * (MU_0/(4*pi)).*dPzdt.*green_kernel;

                %
                %                     interaction_v(time_ind,:) = trapz(d_zprime,trapz(d_xprime, (trapz(d_yprime,(dPhi-electron_velocity.*dA),2) + dPhi_Y0) ,1),3)...
                %                         + trapz(d_xprime,trapz(d_yprime,dPhi_Z0,2),1);

                interaction_v(time_ind,:) = trapz(d_zprime,trapz(d_xprime, (trapz(d_yprime,(-electron_velocity.*dA),2)) ,1),3);
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
                laser , numerical_parameters)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            const = utils.constants_fundamantal();
            EPSILON_0 = const.('EPSILON_0');
            C = const.('C');
            Q_E = const.('Q_E');

            t_c = discretization.t(discretization.t > discretization.t0);
            t_c_subsampled = t_c(1:numerical_parameters.tc_subsampling:end);

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
            gaussian_laser_spot(:,:,abs(zprime) > discretization.z_max,:) = 0;


            l_tc = length(t_c_subsampled);

            n_exc = laser.excited_carriers(material.alpha, material.hnew);
            omega_y = material.calculate_omega_y(n_exc, YPRIME, gaussian_laser_spot);

            t_r = (1/C)*sqrt((x0-XPRIME).^2+(y0-YPRIME).^2+(Z-ZPRIME).^2);

            alpha = material.alpha;
            gamma = material.gamma;
            v_t = material.v_t;
            phase = material.phase;
            gamma_factor = material.gamma_factor;

            begin = '['; last = '] ... Photodember'; arrow = '>';
            body = '='; empty = ' ';
            disp("Calculating potential photodember ....");
            %             omega_y = 0.05 * omega_y;

            parfor time_ind = 1:l_tc

                body_disp = repmat(body , [1,fix(time_ind / length(t_c_subsampled) * 70)]);
                empty_disp = repmat(empty , [1,70 - fix(time_ind / length(t_c_subsampled) * 70)]);
                disp([begin, body_disp, arrow, empty_disp , last]);
                t_prime = t_c_subsampled(time_ind) - t_r;

                g = (Q_E*alpha*v_t^2*...
                    n_exc.*gaussian_laser_spot.*exp(-alpha.*YPRIME)./(omega_y.*(gamma^2+4.*omega_y.^2))).*...
                    (4.*omega_y.*cos(omega_y.*t_prime + phase)+2*gamma.*sin(omega_y.*t_prime + phase)).*...
                    exp(-gamma.*gamma_factor.*t_prime./2);
                g(t_prime<0) = 0;

                dv = green_kernel.*g.*(y0-YPRIME);

                interaction_v(time_ind,:) = (1/(4*pi*EPSILON_0)).* trapz(zprime,trapz(yprime,trapz(xprime,dv,1),2),3);
                %                 Phi(TimeInd,:) = (1/(4*pi*epsilon0)).* trapz(zprime,trapz(yprime,trapz(xprime, GreenKernel.*g.*(y0-Yprime) ,1),2),3);
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


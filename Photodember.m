classdef Photodember
    %PHOTODEMBER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        tc_subsampled
        n_exc
    end

    methods
        function self = Photodember(material,laser,discretization)
            %PHOTODEMBER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function v = interaction_potential(self,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            const = utils.constants_fundamantal();
            HBAR = const.('HBAR');
            Q_E = const.('Q_E');
            EPSILON_0 = const.('EPSILON_0');
            C = const.('C');

            v = zeros(length(self.tc_subsampled),length(self.observation.z));

            tc_sub = self.tc_subsampled;
            green_kernel = self.calculate_green_kernel(self.observation, self.source);
            gaussian_laser_spot = self.laser.gaussian_laser_spot;
            
            l_tc = length(self.tc_subsampled);

            parfor time_ind = 1:l_tc

                t_prime = tc_sub(time_ind) - (1/C)*sqrt((x0-Xprime).^2+(y0-Yprime).^2+(Z-Zprime).^2);
                g = self.name_not_decided(t_prime);
                dv = gaussian_laser_spot.*green_kernel.*g.*(y0-Yprime);
                v(time_ind,:) = (1/(4*pi*EPSILON_0)).* trapz(zprime,trapz(yprime,trapz(xprime,dv,1),2),3);

            end

            %% Return to original t
            v = interp1(self.tc_subsampled,v,tC);
            v = v';
            v = [zeros(length(z),length(t)-length(tC)),v];
            v(isnan(v)) = 0;

        end

        function g = name_not_decided(self, t_prime)

            omega_y = self.material.calculate_omega_y(self, self.material.n_exc, yprime_grid);
            g = (e*self.material.alpha*self.material.v_t^2*...
                self.n_exc.*exp(-self.material.alpha.*Yprime)./(omega_y.*(gamma^2+4.*omega_y.^2))).*...
                (4.*omega_y.*cos(omega_y.*t_prime)+2*gamma.*sin(omega_y.*t_prime)).*...
                exp(-gamma.*t_prime./2);
            g(t_prime<0) = 0;
        end

    end

    methods(Static)

        function green_kernel = calculate_green_kernel(observation, source)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            green_kernel = ((observation.x0-source.xprime).^2 + ...
                (observation.y0-source.yprime).^2+(observation.z-source.zprime).^2).^(-3/2);

        end

    end


end


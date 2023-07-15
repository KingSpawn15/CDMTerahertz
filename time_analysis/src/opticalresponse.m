classdef opticalresponse

    properties 

    end

    methods(Static)

        function nopt = nopt_znte_si(lambda)
            
            lambda = lambda * 1e6;
            nopt = sqrt(4.27 + (3.01*lambda.^2)./(-0.142 + lambda.^2));

        end

        function nopt = nopt_znte(lambda)

            nopt = sqrt(4.27 + (3.01*lambda.^2)./(-0.142 + lambda.^2));

        end

        
        function ngopt = ngopt_znte_si(lambda)
            
            lambda = lambda * 1e6;
            ngopt = -0.5*(lambda*((-6.02*lambda.^3)/(-0.142 + lambda.^2).^2 + ...
                (6.02*lambda)/(-0.142 + lambda.^2)))/sqrt(4.27 + ...
                (3.01*lambda.^2)/(-0.142 + lambda.^2)) + sqrt(4.27 + ...
                (3.01*lambda.^2)/(-0.142 + lambda.^2));
            
        end

        function ngopt = ngopt_znte(lambda)

            ngopt = -0.5*(lambda*((-6.02*lambda.^3)/(-0.142 + lambda.^2).^2 + ...
                (6.02*lambda)/(-0.142 + lambda.^2)))/sqrt(4.27 + ...
                (3.01*lambda.^2)/(-0.142 + lambda.^2)) + sqrt(4.27 + ...
                (3.01*lambda.^2)/(-0.142 + lambda.^2));
            
        end

        function nTHz = nTHz_znte(nu)

            nTHz = sqrt((289.27 - 6*nu.^2)./(29.16 - nu.^2));
            
        end
        
        function nTHz = nTHz_znte_si(nu)
            
            nu = nu * 1e-12;
            nTHz = sqrt((289.27 - 6*nu.^2)./(29.16 - nu.^2));
            
        end
    

        function nTHz = nTHz_inas_drude(nu)
            
            omega = 2 * pi * nu;
            lambda_TO_cm = 217.3;
            lambda_LO_cm = 238.5;
            Ccm = 3e10;
            omega_LO = 2 * pi * Ccm * lambda_LO_cm;
            omega_TO = 2 * pi * Ccm * lambda_TO_cm;
            ninf = 3.5;

            eps_inf = ninf^2;
%             gamma_LO = 2 * pi * Ccm * 2.01;
%             gamma_TO = 2 * pi * Ccm * 8.67;
            gamma_LO = 2 * pi * Ccm * 4.01;
            gamma_TO = 2 * pi * Ccm * 4.01;
            omega_p = 78.1e12;
            gammaD = (1/125e-15);
            background =  eps_inf .* (omega_LO^2 - omega.^2 - 1i * gamma_LO .* omega) ./ ...
                (omega_TO^2 - omega.^2 - 1i * gamma_TO .* omega);

            drude = - 1 .* omega_p^2./(omega.*(omega + 1i*(1/gammaD)));
            nTHz = sqrt(background + drude);
        end


        function nTHz = nTHz_inas(nu)
            
            omega = 2 * pi * nu;
            lambda_TO_cm = 217.3;
            lambda_LO_cm = 238.5;
            Ccm = 3e10;
            omega_LO = 2 * pi * Ccm * lambda_LO_cm;
            omega_TO = 2 * pi * Ccm * lambda_TO_cm;
            ninf = 3.5;

            eps_inf = ninf^2;
            gamma_LO = 2 * pi * Ccm * 2.01;
            gamma_TO = 2 * pi * Ccm * 8.67;

            nTHz = sqrt(eps_inf .* (omega_LO^2 - omega.^2 + 1i * gamma_LO .* omega) ./ ...
                (omega_TO^2 - omega.^2 + 1i * gamma_TO .* omega));
        end

         function nTHz = nTHz_inas_tester(nu)
            
            omega = 2 * pi * nu;
            lambda_TO_cm = 217.3;
            lambda_LO_cm = 238.5;
            Ccm = 3e10;
            omega_LO = 2 * pi * Ccm * lambda_LO_cm;
            omega_TO = 2 * pi * Ccm * lambda_TO_cm;
            ninf = 3.5;
            eps_inf = ninf^2;
            gamma_LO = 2 * pi * Ccm * 2.01;
            gamma_TO = 2 * pi * Ccm * 8.67;
            nTHz = sqrt(eps_inf .* (omega_LO^2 - omega.^2 + 1i * gamma_LO .* omega) ./ ...
                (omega_TO^2 - omega.^2 + 1i * gamma_TO .* omega));
         end

    end
    

end




classdef EELS
    %EELS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        laser
        electron
        discretization
        material
        subsampling
    end
    
    methods
        function self = EELS(electron, laser, discretization , material , subsampling)
            %WAVEFUNCTION Construct an instance of this class
            %   Detailed explanation goes here
            self.laser = laser ;
            self.electron = electron ;
            self.discretization = discretization;
            self.material = material;
            self.subsampling = subsampling;
        end
        
        function interact_v = interaction_v(self, method, interaction_gain_factor,...
                interaction_gain_factor_photodember)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            switch method
                
                case "combination"
                    interact_v = ChargeDynamics.interaction_potential_rectification(self.discretization,...
                        self.material,...
                        self.laser , self.electron, self.subsampling) + ...
                        ChargeDynamics.interaction_potential_photodember(self.discretization, self.material,...
                        self.laser , self.subsampling) * interaction_gain_factor_photodember;
                    
                case "photodember"
                    interact_v = ChargeDynamics.interaction_potential_photodember(self.discretization, self.material,...
                        self.laser , self.subsampling);
                    
                case "rectification"
                    interact_v = ChargeDynamics.interaction_potential_rectification(self.discretization,...
                        self.material,...
                        self.laser , self.electron, self.subsampling);
            end
            
            interact_v = interaction_gain_factor * interact_v;
            interact_v = movmean(movmean(interact_v,3,1),10,2);
            
        end
        
        
        function f_t = calculate_ft(self , interact_v)
            
            t = self.discretization.t;
            omega = self.discretization.omega;
            z = self.discretization.z;
            deltaz = self.discretization.deltaz;
            
            
            %Time domain Fourier transform
            interact_v_fft = fftshift(fft(interact_v,length( t),2),2);
            interact_v_fft = interact_v_fft./max( omega);
            
            % beta(omega)
            [Omega,Z] = meshgrid( omega, z);
            
            interact_v_z = interact_v_fft.*exp(-1i*Omega.*Z/self.electron.electron_velocity);
            beta = sum(interact_v_z,1).* deltaz;
            
            const = utils.constants_fundamantal();
            HBAR = const.('HBAR');
            Q_E = const.('Q_E');
            multiplication_factor = (-1./(1i*HBAR*self.electron.electron_velocity))*Q_E;
            % f(t)
            parfor time_ind = 1:length(t)
                f_t_exp(time_ind) = multiplication_factor.*trapz(omega,2.*real(exp(1i.*omega.*t(time_ind)).*beta));
            end
            
            % Sweep time delay - Psi(E,deltat)
            
            f_t_exp_shift = circshift(f_t_exp,ceil(length(t)/2));
            f_t = exp(-f_t_exp_shift);
            
        end
        
        
        function psi_coherent = calculate_psi_coherent(self, f_t)
            
            deltat = self.discretization.deltat;
            t = self.discretization.t;
            
            delta_t_rep = repmat(deltat',1,length(t));
            t_rep = repmat(t,length(deltat),1);
            
            FT = repmat(f_t,length(deltat),1);
            psi_coherent = FT.*exp(-(t_rep-delta_t_rep).^2/(2*self.electron.electron_time_coherent_sigma.^2));
            
            psi_coherent = fftshift(fft(psi_coherent,length(t),2),2);
            
            psi_coherent = (abs(psi_coherent)).^2;
            psi_coherent = psi_coherent./trapz(self.discretization.energy,psi_coherent,2);
            
        end
        
        
    end
    
    methods(Static)
        function psi_sub = psi_sub_sampled(sub_sample_factor, psi , e_w)
            
            psi_sub = psi(:,1:sub_sample_factor:end);
            psi_sub = psi_sub./trapz(e_w,psi_sub,2);
            
        end
        
        function psi_incoherent = incoherent_convolution(psi, w, t_w, e_w ,...
                w_cut_off_factor)
            
            if nargin < 5
                w_cut_off_factor = 0.01;
            end
            
            psi_sum = zeros(size(psi));
            
            w_cutOff = w_cut_off_factor*max(w(:));
            
            % Example dont use such code. It becomes hard to modify
            % w_cutOff = 0.01*max(w(:));
            
            parfor t_ind = 1:length(t_w)
                for e_ind = 1:length(e_w)
                    
                    if w(t_ind,e_ind) < w_cutOff
                        continue
                    end
                    
                    psi_sum = psi_sum + w(t_ind,e_ind).*...
                        circshift(circshift(psi,-ceil(length(t_w)/2)+t_ind,1),-ceil(length(e_w)/2)+e_ind,2);
                    
                end
            end
            
            psi = psi_sum;
            psi_incoherent = psi./trapz(e_w,psi,2);
            
        end
    end
    
end


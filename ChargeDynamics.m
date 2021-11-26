classdef ChargeDynamics
    %CDEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = ChargeDynamics(inputArg1,inputArg2)
            %CDEM Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function green_kernel = calculate_green_kernel(observation, source)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            green_kernel = ((observation.x0-source.xprime).^2 + ...
                (observation.y0-source.yprime).^2+(observation.z-source.zprime).^2).^(-3/2);

        end
    end
end


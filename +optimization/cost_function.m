function err = cost_function(input_parameters,experimental_parameters)
%COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

err = eels_measurement(input_parameters) - eels_theory(experimental_parameters); 

end


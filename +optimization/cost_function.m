function err = cost_function(x , offset)
%COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    offset = 60;
end

factor_pd_theory = x(1);
factor_or_theory = x(2);
delay = x(3);

list_angle_theory = [0, 90, 30 , 120, 60, 150];
list_angle_measurement = fix(list_angle_theory / 2);

lst = [list_angle_theory;list_angle_measurement];

err = [];

for angle = lst
    
    angle_theory = angle(1);
    angle_measurement = offset + angle(2);
    
    err_i = optimization.error_function(angle_measurement,...
        angle_theory, factor_pd_theory, factor_or_theory , delay);
    
    err = vertcat(err, err_i);
    
end
end


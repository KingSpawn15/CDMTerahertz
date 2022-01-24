function err = cost_function(x , offset)
%COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    offset = 32;
end

factor_pd_theory = x(1);
factor_or_theory = x(2);
delay = x(3);

angle_measurement = offset + 0;
angle_theory = 0;

err1 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory , delay);

angle_measurement = offset + 20/2;
angle_theory = 20;

err2 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory , delay);


angle_measurement = offset + 40 / 2;
angle_theory = 40;

err3 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory, delay);

angle_measurement = offset + 60 / 2;
angle_theory = 60;

err4 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory, delay);

err = [err1; err2; err3;err4];
end


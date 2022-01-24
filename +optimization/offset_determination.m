function err = offset_determination(offset)
%COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

factor_pd_theory = 0.08;
factor_or_theory = 0.14;
delay = 0;


angle_measurement = offset + 0;
angle_theory = 0;

err1 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory , delay);

angle_measurement = offset + 40/2;
angle_theory = 40;

err2 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory , delay);


angle_measurement = offset + 60 / 2;
angle_theory = 60;

err3 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory, delay);

angle_measurement = offset + 90 / 2;
angle_theory = 90;

err4 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory, delay);

err = [err1; err2; err3;err4];


end


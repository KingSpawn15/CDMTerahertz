function err = cost_function(x)
%COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

factor_pd_theory = x(1) * x(3);
factor_or_theory = x(2);

angle_measurement = 178;
angle_theory = 45;

err1 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory);

angle_measurement = 134;
angle_theory = 135;

err2 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory);


angle_measurement = 164;
angle_theory = 75;

err3 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory);

angle_measurement = 118;
angle_theory = 165;

err4 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory);

err = [err1; err2; err3;err4];
end


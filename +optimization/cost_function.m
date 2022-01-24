function err = cost_function(x)
%COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

factor_pd_theory = x(1);
factor_or_theory = x(2);
delay = x(3);

% factor_pd_theory = 0.01;
% factor_or_theory = 0.08;
% delay = x;


angle_measurement = 14;
angle_theory = 0;

err1 = optimization.error_function(angle_measurement,...
    angle_theory, factor_pd_theory, factor_or_theory , delay);

% angle_measurement = 14 + 16 / 2;
% angle_theory = 10;
% 
% err2 = optimization.error_function(angle_measurement,...
%     angle_theory, factor_pd_theory, factor_or_theory , delay);
% 
% 
% angle_measurement = 14 + 44/2;
% angle_theory = 40;
% 
% err3 = optimization.error_function(angle_measurement,...
%     angle_theory, factor_pd_theory, factor_or_theory, delay);
% 
% angle_measurement = 14 + 60 / 2;
% angle_theory = 60;
% 
% err4 = optimization.error_function(angle_measurement,...
%     angle_theory, factor_pd_theory, factor_or_theory, delay);
% 
% err = [err1; err2; err3;err4];

err = err1;
end


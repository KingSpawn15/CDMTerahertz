function err = offset_determination(offset)
%COST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

factor_pd_theory = 0.129;
factor_or_theory = 0.0564;
delay = -25;

% list_angle_theory = [0, 90, 30 , 120, 60, 150];
% list_angle_measurement = fix(list_angle_theory / 2);

list_angle_theory = [100];
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
% angle_measurement = offset + 0;
% angle_theory = 0;
%
% err1 = optimization.error_function(angle_measurement,...
%     angle_theory, factor_pd_theory, factor_or_theory , delay);
%
% angle_measurement = offset + 40/2;
% angle_theory = 40;
%
% err2 = optimization.error_function(angle_measurement,...
%     angle_theory, factor_pd_theory, factor_or_theory , delay);
%
%
% angle_measurement = offset + 60 / 2;
% angle_theory = 60;
%
% err3 = optimization.error_function(angle_measurement,...
%     angle_theory, factor_pd_theory, factor_or_theory, delay);
%
% angle_measurement = offset + 90 / 2;
% angle_theory = 90;
%
% err4 = optimization.error_function(angle_measurement,...
%     angle_theory, factor_pd_theory, factor_or_theory, delay);
%
% err = [err1; err2; err3;err4];


end


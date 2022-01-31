clear all;
% fun = @(x) optimization.cost_function(x);

params.factor_pd_theory = 0.12;
params.factor_or_theory = 0.05;
params.offset = 60;

params.list_angle_theory = 60;
fun_fmincon = @(x) norm(optimization.cost_function_delay(x , params));



lb = [-60];
ub = [60];

intcon = 1;

k = 1;
for angle = 0 : 20 : 180
    params.list_angle_theory = angle;
    x = surrogateopt(fun_fmincon,lb,ub,intcon); 
    delay_struct(k).angle = angle;
    delay_struct(k).delay = x;
    k = k + 1;
end

save('+optimization/optimal_delay.mat','delay_struct');

% close all;
% combination_from_saved_matrices;
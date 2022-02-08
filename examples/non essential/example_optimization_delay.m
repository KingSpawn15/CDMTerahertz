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
for offset = [0, 30, 45, 60, 90]
    for angle = 10 : 10 : 180
        params.list_angle_theory = angle;
        params.offset = offset;
        fun_fmincon = @(x) norm(optimization.cost_function_delay(x , params));
        x = surrogateopt(fun_fmincon,lb,ub,intcon);
        delay_struct(k).offset = offset;
        delay_struct(k).angle = angle;
        delay_struct(k).delay = x;
        k = k + 1;
    end
end

save('+optimization/optimal_delay.mat','delay_struct');

% close all;
% combination_from_saved_matrices;
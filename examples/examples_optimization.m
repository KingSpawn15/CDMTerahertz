clear all;
% fun = @(x) optimization.cost_function(x);
fun_fmincon = @(x) norm(optimization.cost_function(x));

lbarr_pd = [-0.2, -0.1, -0.5];
lbarr_or = [0, 0.6, 0.12];
lbarr_del = [];

lb = [-0.2, 0,  -60];
ub = [0.2, 0.12, 60];

intcon = 3;
x = surrogateopt(fun_fmincon,lb,ub,intcon); 

xarr(1).x = x;



lb = [-0.1, 0,  -35];
ub = [0.1, 0.06, 0];

intcon = 3;
x = surrogateopt(fun_fmincon,lb,ub,intcon); 

xarr(2).x = x;


lb = [-0.1, 0,  0];
ub = [0.1, 0.06, 35];

intcon = 3;
x = surrogateopt(fun_fmincon,lb,ub,intcon); 

xarr(3).x = x;

save('+optimization/last_optimal_xarr.mat','xarr');

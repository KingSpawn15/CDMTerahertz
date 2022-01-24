clear all;
% fun = @(x) optimization.cost_function(x);
fun_fmincon = @(x) norm(optimization.cost_function(x));

% x0 = [0, 0, .50];
lb = [-0.1, 0, -100];
ub = [0.2, 0.6, 100];

% x0 = [50];
% lb = [-75];
% ub = [75];

intcon = 3;
% objconstr = @fun_fmincon;
x = surrogateopt(fun_fmincon,lb,ub,intcon);

% options = optimoptions('lsqnonlin','Display','iter','OutputFcn',@outman);
% [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
% 


clear all;
% fun = @(x) optimization.cost_function(x);
fun_fmincon = @(x) norm(optimization.cost_function(x));

% x0 = [0, 0, .50];
lb = [-0.1, 0, -75];
ub = [0.2, 0.6, 75];

% x0 = [50];
% lb = [-75];
% ub = [75];

intcon = 3;
% objconstr = @fun_fmincon;
x = surrogateopt(fun_fmincon,lb,ub,intcon);

% options = optimoptions('lsqnonlin','Display','iter','OutputFcn',@outman);
% [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
% 
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% 
% 
% nonlcon = @circlecon;
% options = optimoptions('fmincon','Display','iter','OutputFcn',@outman);
% [x,resnorm,residual,exitflag,output] = fmincon(fun_fmincon,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

fun = @(x) optimization.cost_function(x);
x0 = [0.05, 0.1, 0];
lb = [-0.1, 0, -75];
ub = [0.2, 0.6, 75];
options = optimoptions('lsqnonlin','Display','iter');
[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
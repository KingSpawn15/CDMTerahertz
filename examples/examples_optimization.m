fun = @(x) optimization.cost_function(x);
x0 = [0.05, 0.1];
lb = [-0.1, 0];
ub = [0.2, 0.15];
[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub);
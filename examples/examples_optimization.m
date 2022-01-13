fun = @(x) optimization.cost_function(x);
x0 = [0.05, 0.1, 0];
lb = [-0.1, 0, -pi];
ub = [0.2, 0.4, pi];
options = optimoptions('lsqnonlin','Display','iter');
[x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
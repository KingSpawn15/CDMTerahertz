clear all;

lb = [-0.2, -0.2];
ub = [0.2, 0.2];
x0 = [0,0];


delay_vec = -45 : 5 : 0;
offset_vec = [0, 30, 45, 60, 90];
angle_vec = 0 : 180 : 20;
params.list_angle_theory = angle_vec;

for idelay = 1:length(delay_vec)
    for ioffset = 1 : length(offset_vec)
        params.delay = delay_vec(idelay);
        params.offset = offset_vec(ioffset);
        fun = @(x) norm(optimization.cost_function_multi(x , params));
        options = optimoptions('lsqnonlin','Display','iter');
        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options);
        
        index = (idelay-1)*length(offset_vec) + ioffset;
        
        optimation_struct(index).x = x;
        optimation_struct(index).resnorm = resnorm;
        optimation_struct(index).exitflag = exitflag;
        optimation_struct(index).offset = params.offset;
        optimation_struct(index).delay = params.delay;
        
        
    end
end


save('+optimization/optimal_all.mat','optimation_struct');
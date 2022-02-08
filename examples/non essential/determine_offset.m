clearvars;
fun_offset = @(offset) norm(optimization.offset_determination(offset));

err = [];
for offset = mod(0 : 4 : 180 , 180)
    err = vertcat(err,[offset,fun_offset(offset)]);
end

best_offset = sortrows(err , 2);
clearvars;
fun_offset = @(offset) norm(optimization.offset_determination(offset));

err = [];
for offset = mod(-40 : 4 : 40,180)
    err = vertcat(err,[offset,fun_offset(offset)]);
end

best_offset = sortrows(err , 2);
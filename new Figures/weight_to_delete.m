max_eels_sim = [];
max_eels_exp = [];
for inx = 1:11
    max_eels_sim =  horzcat(max_eels_sim, max(eels_comb{inx}));
    max_eels_exp =  horzcat(max_eels_exp, max(com(inx,:)));
end

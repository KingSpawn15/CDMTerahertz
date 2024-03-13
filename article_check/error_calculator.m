function [errs,means] = error_calculator(psi_in,e_w)
%ERROR_CALCULATOR Summary of this function goes here
%   Detailed explanation goes here
errs = []
means = []

for ind = 1 : size(psi_in,1)
    psi_in(psi_in<0) = 0;
    psi = psi_in(ind, :);
    errs = [errs,std(e_w, psi / sum(psi))];
    means = [means,dot(e_w , psi /sum(psi))];
end

end


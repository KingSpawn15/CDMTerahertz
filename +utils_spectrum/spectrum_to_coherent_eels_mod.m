function psi_assemb = spectrum_to_coherent_eels_mod(e_w, t0_vec, loss_spectrum)
    
%     loss_spectrum = interp1(t0_vec.',loss_spectrum.',t_w,'linear',0);
    
    t_w = t0_vec;
    eels_w = loss_spectrum;
    eels_w = repmat(eels_w,[1,length(e_w)]);
    
%     e_w = linspace(-5,5,181);
    e_w_mat = e_w;
    e_w_mat = abs(repmat(e_w_mat,[length(t_w),1]) - eels_w);
    
    [~, b] = min(e_w_mat,[],2);
    
    psi_assemb = zeros(length(t_w), length(e_w));
    row = 1;
    for col = b.'
        psi_assemb(row,col) = 1;
        row = row + 1;
    end

    
end
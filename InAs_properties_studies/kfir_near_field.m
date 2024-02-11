% zzz = -1:0.02:1;
% ttt = -0.1:0.002:0.1;
% [TT, ZZ] = ndgrid(ttt,zzz);

% EOR_k =  kfir_near_field(TT.' * 1e-12 , ZZ.'  * 1e-6, 50e-15, 1e-6);
% close all
% imagesc(TT(:,1), ZZ(1,:), EOR_k)
function EOR_k =  kfir_near_field(TOR, ZOR, TFWHM, x0)

    EOR_k = (-((2*ZOR.^2)./(x0.^2 + ZOR.^2).^2) + 1./(x0.^2 + ...
ZOR.^2))./exp((4*log(2)*TOR.^2)./TFWHM.^2);

% EOR_k = 1./exp((4*log(2)*TOR.*2)./TFWHM.*2);

end

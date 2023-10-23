function [TOR, ZOR, EPDintrap, EOR] = interpolate_field(TOR, ZOR, EOR, TPD, ZPD, EPD)

    EPDintrap = interp2(TPD.', ZPD.', EPD, TOR.', ZOR.', 'linear', 0);

end

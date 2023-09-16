function angle_hwp = angle_hwp_calculator(angle_pol, offset)

    if nargin < 2
        offset = 12;
    end
    
    if mod(offset + fix(angle_pol / 2), 2) == 0 
        angle_hwp = offset + fix(angle_pol / 2);
    else
        angle_hwp = offset + fix(angle_pol / 2) - 1;

    end

end
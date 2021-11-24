function [velocity_msec] = c2msec(velocity_c)

const = utils.constants_fundamantal();
C = const.('C');
velocity_msec = velocity_c * C;

end

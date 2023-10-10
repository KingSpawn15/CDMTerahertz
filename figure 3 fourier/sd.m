function sd()
[mean_eels_measure, deltat, energy, eels_measure, errs] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(90));

end
% 
imagesc(energy, deltat, eels_measure);
hold on;
errorbar(mean_eels_measure, deltat, errs, 'horizontal','r')
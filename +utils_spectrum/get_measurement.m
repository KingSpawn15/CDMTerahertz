function [mean_eels, deltat, energy, eels_measure] = get_measurement(angle)


    [eels_measure, energy,...
        time] = measurement_plots.data_measurement(angle);
    
    deltat = time(time > -1 & time < 1.5);
    eels_measure = eels_measure(time > -1 & time < 1.5,:);
    
    x_c = zeros(1, size(deltat,2));
    
    eels_measure = eels_measure ./ sqrt(sum(eels_measure.^2, 2));
    
    
    for i = 1:size(eels_measure,1)
        [~, ind_max] = max(eels_measure(i,:));
        x_c(i) = ind_max;
    end
    
    mean_eels = energy(floor(x_c));
end

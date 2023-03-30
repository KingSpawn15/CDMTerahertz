[~,~,~] = mkdir('results/measurement_plots_reconfig');
for angle = [0:2:180]
    close all
    measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
    
    offset = 100;
    theta = num2str(mod(fix(2*(-offset + angle)),360));
    
    str = ['results/measurement_plots_reconfig/','eels_angle=',num2str(theta)];
    savefig(gcf,strcat(str,'.fig'));
    exportgraphics(gcf,strcat(str,'.png'),'Resolution',300);
end


[~,~,~] = mkdir('results/measurement_plots');
for angle = [0:2:180]
    close all
    measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);
    str = ['results/measurement_plots/','eels_angle=',num2str(angle)];
    savefig(gcf,strcat(str,'.fig'));
    exportgraphics(gcf,strcat(str,'.png'),'Resolution',300);
end


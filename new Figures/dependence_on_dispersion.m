clear all
close all

lambda = 800e-9;
tau = 30e-15;
z = -1e-6;
d  = .5e-3;

materials = opticalresponse;
refractive_index_data = read_refractive_index('refractive_index_data/InAs.txt');
nopt = @(lambda) interpolate_refractive_index(refractive_index_data, lambda * 1e9);
delta_lambda = 0.1*1e-9;
ngopt = @(lambda) nopt(lambda) - (lambda)*(nopt(lambda + delta_lambda) - ...
    nopt(lambda))/(delta_lambda);
% ngopt = @(lambda) 1
% nopt = @(lambda) 3.5 + 0.0000001i


np = 10001;
omega_max = 2 * pi * 12e12;

i = 1;

time_ps_c = {};
ethz_t_c = {};
for factor = 1
    disp(i)
%     nTHz = @(omega) materials.nTHz_inas_drude_modified(omega / (2 * pi) , (78.1/14.7));
    nTHz = @(omega) 1;
    [time_ps, ethz_t, omg, eomg] = electric_field_time(lambda, tau, z, d, ...
        ngopt, nTHz, nopt, ...
        np, omega_max);
    ethz_t_re = real(ethz_t);
    time_ps_c{i} = time_ps;
    ethz_t_c{i} = ethz_t_re;
    i = i + 1;
end

%%
close all
setdir = 'new Figures/results/';
FontName = 'ariel';
FontSize = 15;
intensity = @(t,maxI) maxI * exp(-t.^2 / (2 * (tau * 1e12).^2));
maxI = 1;
a = area(time_ps - 0.1,intensity(time_ps,maxI));
a.FaceAlpha = 0.5;
a.LineStyle = "none";
a.FaceColor = '#FF8C00'

% hold on
% plot(cell2mat(time_ps_c)-0.1, cell2mat(ethz_t_c) ./ max(cell2mat(ethz_t_c)), ...
%     'Color','#00008B','LineWidth',1); 


set(gca,'FontSize',FontSize);
set(gca,'TickLength',[.02,.02]);
yticklabels([])
xlim([-.3,1.5]);
ylim([-0.5,1.5])
% 
pbaspect([1 1 1])
set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
exportgraphics(gcf, [setdir, 'section_EOR_no_dispersion.png'],'resolution', 300);

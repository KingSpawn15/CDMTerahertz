clear all
close all

fields = optimal_parameters();
EOR = fields.EOR;
EPD = fields.EPD;
e_w = fields.e_w ;
T = fields.T ;
Z = fields.Z ;
eels_photodember = fields.eels_obj;
%% Get measurement data

angle = sort([0, 30, 50, 70, 90, 110, 130, 150, 180]);

% [~ ,~, ~, eels_measure_0] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(0));
for ii = 1:length(angle)
[~, deltat, energy, eels_measure{ii}] = utils_spectrum.get_measurement(utils_spectrum.angle_hwp_calculator(angle(ii)));
end


%%

for ii = 1:length(angle)
    disp(angle(ii))
%     psi_incoherent_comb{ii} = generate_incoherent_spectrum_for_angle(kfactor * EPDintrap, EOR, TOR, ZOR, eels, w, t_w, e_w, angle(ii));
    [e_w, t_w, psi_incoherent_comb{ii}] = utils_spectrum.generate_incoherent_spectrum_for_angle(EPD, EOR, T, Z, angle(ii), eels_photodember, e_w);
end

%%
close all
close all
setdir = 'new Figures/results/';
figure;
image_name = 'tiles_multiple';
FontName = 'ariel';
FontSize = 10;
ttt = tiledlayout(2,length(angle),"TileSpacing","compact");
ttt.Padding = "compact";

for ii = 1:length(angle)
plot_tile(energy, deltat, eels_measure{ii});

set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[],'','',FontSize,[0 0 0, 0]);
end

end

for ii = 1:length(angle)
plot_tile(e_w,t_w, psi_incoherent_comb{ii});
set_axis_properties(gca,FontSize,FontName,0.01,[],[],'','',FontSize,[0 0 0, 0]);
if ii ==1
set_axis_properties(gca,FontSize,FontName,0.01,[-1:0.5:1.5],[-4:2:4],'','',FontSize,[0 0 0, 0]);
end
end


set(gcf,'Position',[200,200,200 + 120 * length(angle), 200 +  160]);


exportgraphics(gcf, [setdir,image_name,'.png'], 'Resolution',300);
%%
function plot_tile(x, y, z)
    FontSize = 10;
    FontName = 'ariel';
    nexttile
    imagesc(x, y, z);
    ylim([-1,1.5]);
    xlim([-5,5]);
    colormap jet
    axis square
    set_axis_properties(gca,FontSize,FontName,1,[],[],'','',FontSize,[0.3 0.3 0.3]);
end

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end


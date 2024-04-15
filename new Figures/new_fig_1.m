clear all
close all

fields = optimal_parameters();
EOR = fields.EOR;
EPD = fields.EPD;
e_w = fields.e_w ;
T = fields.T ;
Z = fields.Z ;
eels_photodember = fields.eels_obj;
%%
close all;
FontName = 'ariel';
FontSize = 15;
clim  = max(abs(EOR(:)));
clim_pd  = max(abs(EPD(:)));
setdir = 'new Figures/results/';
create_figure_electricfield(T+.1, Z, EOR, clim, setdir, 'field_rectification.png', FontSize);
create_figure_electricfield(T+.1, Z, movmean(movmean(EPD,5,2),2,1), clim_pd, setdir, 'field_photodember.png', FontSize);

%%
close all
figure;
plot(T(:,1), EOR(401,:),'Color','#00008B','LineWidth',1)
xlim([-0.3,1.5]); %ylim([-3,3]*1e6)
tau = 30e-15 * 1e12;
intensity = @(t,maxI) maxI * exp(-t.^2 / (2 * tau.^2)); 
maxI = max(EOR(401,:));
hold on
a = area(T(:,1) - 0.1,intensity(T(:,1),maxI));
a.FaceAlpha = 0.4;
a.LineStyle = "none"
pbaspect([1 1 1])
set(gca,'FontSize',FontSize);
set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
exportgraphics(gcf, [setdir, 'section_EOR.png'],'resolution', 300);


figure;
plot(T(:,1), EPD(412,:),'Color','#FFA500','LineWidth',1)
xlim([-0.3,1.5]); %ylim([-5,5]*1e4)
tau = 30e-15 * 1e12;
intensity = @(t,maxI) maxI * exp(-t.^2 / (2 * tau.^2)); 
maxI = 2644720;
pbaspect([1 1 1])
set(gca,'FontSize',FontSize);
set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
exportgraphics(gcf, [setdir, 'section_EPD.png'],'resolution', 300);
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


function ll = create_line(deltat)
    c = 3*10^(8 - 12 +6);
    ve = 0.7*c;
    z = -100 : 100;
    t = deltat + z / ve;
    ll = line(t,z,'LineWidth',2, 'Color',[0,0.7,0]);
    ll.Color = [ll.Color 0.4];
end

function create_figure_electricfield(T, Z, E, clim, setdir, filename, FontSize)
    figure;
    imagesc(T(:,1), Z(1,:), E, 'Interpolation', 'bilinear',[-clim, clim]);
    set(gca,'FontSize',FontSize);
    xlim([-.3,1.5]);
    ylim([-60,60]);
    yticks(-60:20:60)
    xticks(-.3:.3:1.5)
    colormap(utils.redblue);
    pbaspect([1 1 1])
    set(gcf,'position', [200 , 200 , 200 + 200, 200 + 120]);
    colorbar;
    l0 = create_line(0); l1 = create_line(0.3); l2 = create_line(0.7); 

    set(groot,'defaultAxesXTickLabelRotationMode','manual');
    exportgraphics(gcf, [setdir, filename],'resolution', 300);

    
end

function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
%     axes tight;
end
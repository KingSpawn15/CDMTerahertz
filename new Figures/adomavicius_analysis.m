xdata1 = [10, 20, 40, 60, 80, 115, 160, 200];
xdata2 = [10, 20, 40, 55, 80, 120, 160, 200];
ydata1 = [40, 85, 135, 190, 225, 300, 370, 400];
ydata2 = [30, 60, 70, 85, 100, 115, 130, 130];

power = @(param,xdata) param(1).*(xdata.^param(2));

saturate = @(param,xdata) param(1).*(xdata)./(xdata + abs(param(2)));

% Initial guess for the parameters [A, Gamma, x0]
param0 = [200, 30];
param02 = [10, 43.7523];
% Your xdata and ydata here
% xdata = pulse_energy_list(:);
% ydata = eor_max;

% Use lsqcurvefit to fit the Lorentzian to your data
options = optimset('Display','off');
[param1,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(saturate, param0, xdata1, ydata1, [], [], options);
[param2,resnorm2,residual2,exitflag2,output2,lambda2,J2] = lsqcurvefit(saturate, param02, xdata2, ydata2, [], [], options);


N = length(xdata1) - length(param1);
stdErrors = full(sqrt(diag(inv(J'*J)*resnorm/N)));

N2 = length(xdata2) - length(param2);
stdErrors2 = full(sqrt(diag(inv(J2'*J2)*resnorm/N2)));


% close all
% figure;
% ttt = tiledlayout(1,1,"TileSpacing","compact");
% ttt.Padding = "loose";
% nexttile
% plot(pulse_energy_list, eor_max)
% nexttile
% plot(pulse_energy_list, eor_max)
% hold on;
%% 
f1 = @(x)saturate(param1,x); hold on;
f2 = @(x)saturate(param2,x);
xx = 0:5:max(xdata1);
fit1 = f1(xx);
fit2 = f2(xx);

close all
t = tiledlayout(1,1);
ax1 = axes(t);
ax2 = axes(t);
ps = plot(ax2, ...
    xdata1 * 0.93 ,ydata1, ...
    xdata2 * 0.93, ydata2,...
    'LineStyle','none','MarkerSize',10); hold on;

set(ps, ...
    {'Marker'},{'s';'o'}, ...
    {'MarkerFaceColor'},{'#5928ED';'#5BA300'}, ...
    {'MarkerEdgeColor'},{'none';'none'}...
    )

ps2 = plot(ax1,xx, fit1, xx, fit2, 'LineWidth',1); 
set(ps2, ...
    {'Color'},{'#5928ED';'#5BA300'});
xlim(ax1,[0,225]);
ylim(ax1,[0,450]);
xlim(ax2,[0,225]);
ylim(ax2,[0,450]);

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.XColor = '#b30000';
ax2.YTick = [];
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

set_axis_properties(ax1,FontSize+4,FontName,1,[0:100:400],0:50:200,'','',FontSize,[0.3 0.3 0.3]);
set_axis_properties(ax2,FontSize+4,FontName,1,[],[0:50:200],'','',FontSize,[0.3 0.3 0.3]);



% axis off;
set(gca, 'Color', 'None');
set(gcf,'Position',[200,50,200 + 600,200 + 300]);
exportgraphics(gcf,'new Figures/results/adom_compare.PNG','BackgroundColor','none')

%%
function set_axis_properties(ax,FontSize,FontName,LineWidth,YTick,XTick,ylabel_str,xlabel_str,label_FontSize,label_Color)
    ax.FontSize = FontSize;
    ax.FontName = FontName;
%     ax.LineWidth = LineWidth;
    ax.YTick = YTick;

    ax.XTick = XTick;

    ylabel(ylabel_str,'Color',label_Color,'FontSize',label_FontSize);
    xlabel(xlabel_str,'Color',label_Color,'FontSize',label_FontSize);
end
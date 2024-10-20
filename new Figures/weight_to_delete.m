max_eels_sim = [];
max_eels_exp = [];
for inx = 1:11
    max_eels_sim =  horzcat(max_eels_sim, max(eels_comb{inx}));
    max_eels_exp =  horzcat(max_eels_exp, max(com(inx,:)));
end


%%
% Data for the first two plots
% x = linspace(0, 10, 100);
close all
y1 = sin(x);
y2 = cos(x);

% Data for the fplot
f = @(x) exp(-x).*sin(2*pi*x);

% Data for the last plot
y3 = x.^2;

% Create a tiled layout with 1 tile
T = tiledlayout(1, 1);

% First axis: y1 and y2 with x-axis on top
ax1 = axes(T);
plot(ax1, pulse_energy_list * 159, ecomb_max, 'ro');
hold on;
plot(ax1, pulse_energy_list * 159, eor_max, 'bo');
plot(ax1, pulse_energy_list * 159, epd_max * 10, 'go');
set(ax1, 'XAxisLocation', 'top'); % Move x-axis to the top
xlabel(ax1, 'X-axis on top');
ylabel(ax1, 'Y1 and Y2 axis');
title(ax1, 'Two plots with same x-axis');

% Second axis: fplot with different x-axis and no y-ticks
ax2 = axes(T);
% fplot(ax2, f, [0 2*pi], 'g');
fplot(ax2,@(x)saturate(param,x), [min(xdata), max(xdata)],'color','#b30000',LineWidth=2);
set(ax2, 'YAxisLocation', 'left', 'YTick', [], 'Color', 'none'); % No y-ticks, transparent background
% ax2.Position = ax1.Position; % Align with the first axis
ax2.XAxisLocation = 'bottom'; % Ensure the x-axis is at the bottom
xlabel(ax2, 'Different X-axis');

% % Third axis: y3 with y-axis on the right
ax3 = axes(T);
plot(ax3, pulse_energy_list, max_exp, 's');
yyaxis(ax3, 'right'); % Y-axis on the right
set(ax3, 'XAxisLocation', 'bottom', 'XTick', [], 'Color', 'none'); 
% ax3.Position = ax1.Position; % Align with the first axis
ax3.XAxis.Visible = 'off'; % Hide the x-axis since it overlaps
ylabel(ax3, 'Y3 axis (right)');

% % Ensure all axes are linked
% linkaxes([ax1, ax2], 'x'); % Link x-axes

ylim(ax2, [0,4.5e6]);
ylim(ax1, [0,4.5e6]);

set_axis_properties(ax1,FontSize+4,FontName,1,1e6*[0:0.5:6],0:5:30,'','',FontSize,[0.3 0.3 0.3])
set_axis_properties(ax2,FontSize+4,FontName,1,1e6*[0:0.5:6],[0:5:30] * 200,'','',FontSize,[0.3 0.3 0.3])
%%
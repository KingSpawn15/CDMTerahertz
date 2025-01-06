% Specifics pertaining to the measurements
clear all
close all
%% Get measurement data
[com,deltat,energy,eels_measure,errs] = getmeasurement_power_2();
[com_2,deltat_2,energy_2,eels_measure_2,errs_2] = getmeasurement_power();
% get the spectrum for all powers
pulse_energy_list = [0.1000    0.1700    0.3000    0.5200    1.0000    1.7300    3.0000    5.1900   10.0000 17.3 30];

%% plot experimental max loss for all powers
max_loss = zeros(11,1);
ind = 1:11
for ind = 1:11
    max_loss(ind) = max(com_2(ind,2:1:end));

end

saturate = @(param,xdata) param(1).*(xdata)./(xdata + abs(param(2)));

% Initial guess for the parameters [A, Gamma, x0]
param0 = [1, 5];

% Your xdata and ydata here
xdata = pulse_energy_list(1:end).';
ydata = max_loss(1:end);

% Use lsqcurvefit to fit the Lorentzian to your data
options = optimset('Display','off');
[param,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(saturate, param0, xdata, ydata, [], [], options);


N = length(xdata) - length(param);
stdErrors = full(sqrt(diag(inv(J'*J)*resnorm/N)));

plot(pulse_energy_list, max_loss)

close all
figure;
plot(pulse_energy_list, max_loss,'LineStyle','none', ...
    'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
hold on;
fplot(@(x)saturate(param,x), [min(xdata), max(xdata)],LineWidth=2);
fplot(@(x) param(1) .*x ./ param(2), [min(xdata), max(xdata)],'--',LineWidth=2);
xlim([-5, 35]);
ylim([0, param(1) * 1.5]);
%%
% error analysis figure
close all
figure;
set(groot,'defaultAxesXTickLabelRotationMode','manual')
FontName = 'ariel';
FontSize = 14;
ttt = tiledlayout(1,2,"TileSpacing","compact");
% ax1 = axes(ttt);
% ax1.Layout.Tile = 5;
ttt.Padding = "loose";
nexttile
imagesc(energy, deltat- 0.3, eels_measure{10});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,-1:0.5:1.5,-4:2:4,'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
hold on
errorbar(com_2(10,2:1:end),deltat(2:1:end) - .3,errs(10,2:1:end)/2,'horizontal', ...
    'Color',[1 1 1] * .7, ...
    'LineWidth',1.5,'LineStyle','none');

nexttile
imagesc(energy, deltat- 0.3, eels_measure{10});
    ylim([-1,1.5]);
    xlim([-5,5]);
    set_axis_properties(gca,FontSize,FontName,1,[],-4:2:4,'','',FontSize,[0.3 0.3 0.3])
    colormap jet
    axis square
hold on
errorbar(com(10,2:1:end),deltat(2:1:end) - .3,errs(10,2:1:end)/2,'horizontal', ...
    'Color',[1 1 1] * .7, ...
    'LineWidth',1.5,'LineStyle','none');
set(gcf,'Position',[200,200,200 + 400,200 + 200]);
exportgraphics(gcf, 'article_check/results/error_analysis.png', 'Resolution',300);

max_exp = []
for i = 1:11
    max_exp = [max_exp, max(abs(com_2(i,:)))]
end


%%
function [errorx, errory] = errorband(mean, error, deltat)
    errorx = [mean + error/2, flip(mean - error/2)];
    errory = [deltat, flip(deltat)];
end


function [com,deltat, energy, eels_measure,errs] = getmeasurement_power()
    load('saved_matrices\PulseEnergy.mat')
    deltat = Time;
    energy = EnergyCrop;

    
    x_c = zeros(11, size(Time,2));
    errs = ones(11, size(Time,2));

    for ii = 1:11
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        eels_measure{ii} = eels;
         
        mass_matrix = eels;
        % Initialize a vector to store the x coordinates of the center of mass of each row
%         errs = zeros(1, size(deltat,2));
        for i = 1:size(eels,1)
%             row_mass = sum(mass_matrix(i,:));
            for j = 1:size(eels,2)
                [~, ind_max] = max(mass_matrix(i,:));
                x_c(ii,i) = ind_max;
            end
            x_c(ii,i) = x_c(ii,i);

            errs(ii,i) = energy(find(mass_matrix(i,:) > 0.7 * max(mass_matrix(i,:)),1,'last')) - ...
            energy(find(mass_matrix(i,:) >  0.7 * max(mass_matrix(i,:)),1,'first'));

        end
    end

    com = energy(floor(x_c));
end


function [com,errs] = get_errors_sim(eels_sim, e_w, t_w)
        
    com = zeros(11, size(t_w,1));
    errs = ones(11, size(t_w,1));
    
    for ii = 1:11
        eels = eels_sim{ii};
        eels = eels ./ sqrt(sum(eels.^2, 2));
        [errs_i,means_i] = error_calculator(eels,e_w);

        com(ii,:) = means_i;
        errs(ii,:) = errs_i;
    
    end

end




function [com,deltat, energy, eels_measure,errs] = getmeasurement_power_2()
    load('saved_matrices\PulseEnergy.mat')
    
    deltat = Time;
    energy = EnergyCrop;
    
    com = zeros(11, size(Time,2));
    errs = ones(11, size(Time,2));

    for ii = 1:11
        
        eels = squeeze(cell2mat(DataSetCropAll(ii)))';
        eels = eels ./ sqrt(sum(eels.^2, 2));
        [errs_i,means_i] = error_calculator(eels,energy);

        eels_measure{ii} = eels;
        com(ii,:) = means_i;
        errs(ii,:) = errs_i;
    
    end
% 
%     [errs,means] = error_calculator(psi_in,e_w)
%     
%     x_c = zeros(11, size(Time,2));
%     errs = ones(11, size(Time,2));
% 
%     for ii = 1:11
%         eels = squeeze(cell2mat(DataSetCropAll(ii)))';
%         eels = eels ./ sqrt(sum(eels.^2, 2));
%         eels_measure{ii} = eels;
%          
%         
% 
%         end
%     end
% 
%     com = energy(floor(x_c));
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

function [errs,means] = error_calculator(psi_in,e_w)
%ERROR_CALCULATOR Summary of this function goes here
%   Detailed explanation goes here
errs = []
means = []

for ind = 1 : size(psi_in,1)
    psi_in(psi_in<0) = 0;
    psi = psi_in(ind, :);
    errs = [errs,std(e_w, psi / sum(psi))];
    means = [means,dot(e_w , psi /sum(psi))];
end

end


clear all;
close all;
dir = 'results/fourier_analysis/';
[~,~,~] = mkdir(dir);
ii = 1;


theta_chart = [];
for angle = 0 : 2: 178
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    offset = 100;
    theta_chart = vertcat(theta_chart,[mod(fix(2*(-offset + angle)),360),angle]);

end
theta_chart = sortrows(theta_chart);

ii = 1;
for angle = [0: 2: 0]
    %     angle = 10;
    %     [EnergyCrop,time,dataset] = measurement_plots.plot_measurement('saved_matrices/rectification.mat', angle);

    [dataset, EnergyCrop,...
        time] = measurement_plots.data_measurement(angle);

    offset = 100;
    theta = mod(fix(2*(-offset + angle)),360);

    if ~((mod(theta,10) == 0) || (mod(theta+2,10) == 0))
        continue
    end

    theta_str = num2str(theta);

    tc = time(time > -1 & time < 1.5);
    psi = dataset(time > -1 & time < 1.5,:);
    emat = repmat(EnergyCrop,[length(tc),1]);
    psi_e = psi.*emat;
    psi_exp_e = trapz(EnergyCrop,psi_e,2)./trapz(EnergyCrop,psi,2);
    psi_omega = fft(fftshift(psi_exp_e));

    dt = (tc(2) - tc(1))*1e-12;
    fs = 1/dt;
    lsignal = length(tc);
    Y = fft(psi_exp_e);
    P2 = abs(Y/lsignal);
    P1 = P2(1:lsignal/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(lsignal/2))/lsignal;

    str_str = strcat('angle',num2str(theta));
    data.x = f*1e-12;
    data.y = P1;
    fft_struct.(str_str) = data;
    ii = ii + 1;
%%
    close all;
    figure;
    set(gcf,'Position',[100,50,600,450]);

    plot(f*1e-12,P1,'-o','LineWidth',3,'Color','#6685E1');
    xlim([0,5]);
    ylim([0,1.2]);
    set(gca,'FontSize',14);
    xlabel('Frequency (THz)','FontSize',18,'FontName','Times New Roman');
    ylabel('Power','FontSize',18,'FontName','Times New Roman');

    hold on;
    axes('pos',[.45 .4 .45 .5]);
    imagesc(EnergyCrop,tc,psi);
    annotation('textbox',[0.45, 0.8, 0.1, 0.1],'String',['\theta = ',theta_str,'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',18);
    annotation('textbox',[0.45, 0.4, 0.1, 0.1],'String',['\Lambda = ',num2str(angle),'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',8);
    colormap('bone');
    colorbar;
    set(gca,'FontSize',10);
    xlabel('Energy shift [eV]','FontSize',14,'FontName','Times New Roman');
    ylabel('Time delay [ps]','FontSize',14,'FontName','Times New Roman');
    yticks([-1 : 0.2 : 1.5]);

    xticks([-4:2:4]);
    ylim([-1,1.5]);
    xlim([-5,5]);

%%
    str = [dir,...
        'single_theta=',num2str(theta),...
        ];
    exportgraphics(gcf, strcat(str,'.png'),'resolution' , 400);


end


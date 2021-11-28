clc; close all; clear all;

%% Parameters
TotalRepetitions = 3;
Time = -0.75:0.03:3;%[ps]
PulseEnergyNum = 11;

PathMain = 'measurement_data\111 Power scans\';
MainFolderList = dir(PathMain);

for i = 1:PulseEnergyNum
    PulseList{i} = MainFolderList(i+2).name(1:end-6);
    PulseListNum(i) = str2double(MainFolderList(i+2).name(1:end-8));
end
[~,Indx] = sort(PulseListNum);
PulseList = PulseList(Indx);
PulseListNum = PulseListNum(Indx);
PulseListLength = length(PulseList);

%% Read .dm4 files
for PulseEnergyInd = 1:PulseListLength
    
    CurrentFolderList = dir([PathMain,PulseList{PulseEnergyInd},' 58deg']);
    CurrentFolderHeading = CurrentFolderList(3).name(1:14);

    for RepetitionInd = 1:TotalRepetitions            
        
        FilePath = [PathMain,PulseList{PulseEnergyInd},' 58deg','\',CurrentFolderHeading,'_00',num2str(RepetitionInd)];
        saveFlag = false;
        FileName = nan;
                              
        RepNumStr = num2str(RepetitionInd);
        [DataSet,Energy,Units] = Readdm4_111PowerScans(FilePath,RepNumStr);
        Position = 0;% Only one position at the sample face
        
        %Crop around ZLP area
        [~,EnergyZeroInd] = max(DataSet(:,:,1),[],2);
        EnergyZero = Energy(EnergyZeroInd);

        EnergyDelta = 5;%[eV]
        EnergyMax = EnergyZero + EnergyDelta;
        EnergyMin = EnergyZero - EnergyDelta;
        
        for PosInd = 1:length(Position)
            EnergyCrop(PosInd,:) = Energy(Energy >= EnergyMin(PosInd) & Energy <= EnergyMax(PosInd))...
                                    - EnergyZero(PosInd);
            for TimeInd = 1:length(Time)
                dataset = DataSet(PosInd,:,TimeInd);
                DataSetCrop(PosInd,:,TimeInd) = dataset(Energy >= EnergyMin(PosInd) & Energy <= EnergyMax(PosInd));
            end
        end

        EnergyCropAll{PulseEnergyInd,RepetitionInd} = EnergyCrop;    
        DataSetCropAll{PulseEnergyInd,RepetitionInd} = DataSetCrop;

    end
end

%% Analysis
for PulseEnergyInd = 1:PulseListLength
    for RepetitionInd = 1:TotalRepetitions          
        for PosInd = 1:length(Position)

            EnergyCrop = squeeze(EnergyCropAll{PulseEnergyInd,RepetitionInd}(PosInd,:))';        
            DataSetCrop = squeeze(DataSetCropAll{PulseEnergyInd,RepetitionInd}(PosInd,:,:))';

            CoG = trapz(EnergyCrop,EnergyCrop'.*DataSetCrop,2)./trapz(EnergyCrop,DataSetCrop,2);
            
            %% Denoise
            CoG = movmean(CoG,3);        
            CoG = CoG - mean(CoG);
            
            %%
            CoGAll{PulseEnergyInd,RepetitionInd}(PosInd,:) = CoG;
            
        end
    end
end

%% Plot deltat-position maps & cross-sections
close all;

for PulseEnergyInd = 1:PulseListLength

    for RepetitionInd = 1:TotalRepetitions          

        figure(PulseEnergyInd)
        subplot(1,TotalRepetitions,RepetitionInd)
        plot(Time,CoGAll{PulseEnergyInd,RepetitionInd},'linewidth',2)
        xlabel('Time delay [ps]')
        ylabel('Energy shift [eV]')
        xlim([min(Time) max(Time)])        
        axis square
        ax = gca;
        ax.FontSize = 18;
        ax.LineWidth = 1;
        box on
        
        CoGAllRep(RepetitionInd,:,:) = CoGAll{PulseEnergyInd,RepetitionInd};
                
        PV(PulseEnergyInd,RepetitionInd) = max(CoGAllRep(RepetitionInd,:,:)) - min(CoGAllRep(RepetitionInd,:,:));

    end
    
    CoGAllRepMean = squeeze(mean(CoGAllRep,1));
    CoGAllRepSTD = squeeze(std(CoGAllRep,0,1));

    figure(PulseEnergyInd*1000)
    plot(Time,CoGAllRepMean,'linewidth',2)
    xlabel('Time delay [ps]')
    ylabel('Energy shift [eV]')
    xlim([min(Time) max(Time)])
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    box on

    [~,MaxTimeInd] = max(CoGAllRepMean);
    
    time = Time - Time(MaxTimeInd);
    if PulseEnergyInd == 1
        time = Time;
    end
    
    colormap jet
    cmap = colormap;
    cmap = cmap(1:round(length(cmap)/PulseListLength):end,:);
%     cmap = flip(cmap,1);

    figure(100)
    hold on
    errorbar(time,CoGAllRepMean + 1.5*(PulseEnergyInd-1),CoGAllRepSTD,'linewidth',2,'color',cmap(PulseEnergyInd,:))
    hold off
    xlabel('Time delay [ps]')
    ylabel('Energy shift [eV]')
    xlim([min(Time) 2.7])
    ylim([-1 18])
    legend(flip(PulseList))
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    box on
    
end

for PulseEnergyInd = 1:PulseListLength
    figure(100)
    hold on
    plot([-0.75 3],1.5*(PulseEnergyInd-1)*[1 1],'--k','linewidth',1)
    hold off
end

figure(200)
hold on
errorbar(PulseListNum,mean(PV,2),std(PV,0,2),'xk','linewidth',2)
hold off
xlabel('Laser pulse energy [nJ]')
ylabel('P-V energy shift [eV]')
xlim([0 30])
ylim([0 4])
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
box on

%% Plot time scans
close all;

for PulseEnergyInd = 1:PulseListLength
    for RepetitionInd = 1:TotalRepetitions          
        
        DataSetAllRep(RepetitionInd,:,:) = squeeze(DataSetCropAll{PulseEnergyInd,RepetitionInd}(1,:,:));
        Energy = squeeze(EnergyCropAll{PulseEnergyInd,RepetitionInd}(1,:));
        
    end
    
    dataset = squeeze(mean(DataSetAllRep,1));
    dataset = movmean(dataset,3,2);
    dataset = dataset./trapz(Energy,dataset',2)';

    figure(PulseEnergyInd)
    imagesc(Energy,Time,dataset')
    xlabel('Energy shift [eV]')
    ylabel('Time delay [ps]')
    xlim([min(Energy) max(Energy)])
    ylim([min(Time) max(Time)])
    colormap jet
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    box on
    caxis([0 max(dataset(:))])

end

%% Plot Fourier analysis
t = [Time,3.03:0.03:10];

dt = (t(2)-t(1))*1e-12;
fs = 1./dt;
L = 2*length(t);
df = fs/L;%[Hz]
t = (-L/2:L/2)*dt;%[s]
f = (-L/2:L/2)*df;%[Hz]

% s = sin(2*pi*5e12.*t);
% figure
% plot(s)

for PulseEnergyInd = 1:PulseListLength
    Fs = zeros(TotalRepetitions,length(t));
    for RepetitionInd = 1:TotalRepetitions          
       
        s = CoGAll{PulseEnergyInd,RepetitionInd}(1,:);
        s = s./max(s);
        s = [zeros(1,(length(t)+1)/2),s,zeros(1,(length(t)-1)/2-length(s))];

        Fs(RepetitionInd,:) = fftshift(fft(s));
        
    end
    
    Fs = abs(Fs)./max(abs(Fs),[],2);

    FsMean = mean(Fs,1);
    FsSTD = std(Fs,0,1);
    
    figure(300)
    hold on
    errorbar(f*1e-12,FsMean + (PulseEnergyInd-1),FsSTD,'linewidth',2,'color',cmap(PulseEnergyInd,:))
    hold off
    xlabel('Frequency [THz]')
    ylabel('Power spectral density (Norm.)')
    xlim([0 5])
    ylim([-0.5 11.5])
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    box on
end

for PulseEnergyInd = 1:PulseListLength
    figure(300)
    hold on
    plot([0 5],(PulseEnergyInd-1)*[1 1],'--k','linewidth',1)
    hold off
end
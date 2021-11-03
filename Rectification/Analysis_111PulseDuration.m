clc; close all; clear all;

%% Parameters
TotalRepetitions = 3;
Time = -0.75:0.03:3;%[ps]
PulseDurationNum = 11;%

PathMain = 'C:\Users\yanna\Dropbox\PhD\Research\Photo-Dember\Experiment\Session 2\111 Pulse duration\';
MainFolderList = dir(PathMain);

for i = 1:PulseDurationNum
    HeaderInd = find(MainFolderList(i+2).name=='_');
    PulseListFullName{i} = MainFolderList(i+2).name;
    PulseList{i} = MainFolderList(i+2).name(HeaderInd+1:end);
    PulseListNum(i) = str2double(MainFolderList(i+2).name(HeaderInd+1:end-2));
end
[~,Indx] = sort(PulseListNum);
PulseList = PulseList(Indx);
PulseListFullName = PulseListFullName(Indx);
PulseListNum = PulseListNum(Indx);
PulseListLength = length(PulseList);

%% Read .dm4 files
for PulseDurationInd = 1:PulseListLength
    
    CurrentFolderList = dir([PathMain,PulseListFullName{PulseDurationInd}]);
    CurrentFolderHeading = CurrentFolderList(3).name(1:14);

    for RepetitionInd = 1:TotalRepetitions            
        
        FilePath = [PathMain,PulseListFullName{PulseDurationInd},'\',CurrentFolderHeading,'_00',num2str(RepetitionInd)];
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

        EnergyCropAll{PulseDurationInd,RepetitionInd} = EnergyCrop;    
        DataSetCropAll{PulseDurationInd,RepetitionInd} = DataSetCrop;

    end
end

%% Analysis
for PulseDurationInd = 1:PulseListLength
    for RepetitionInd = 1:TotalRepetitions          
        for PosInd = 1:length(Position)

            EnergyCrop = squeeze(EnergyCropAll{PulseDurationInd,RepetitionInd}(PosInd,:))';        
            DataSetCrop = squeeze(DataSetCropAll{PulseDurationInd,RepetitionInd}(PosInd,:,:))';

            CoG = trapz(EnergyCrop,EnergyCrop'.*DataSetCrop,2)./trapz(EnergyCrop,DataSetCrop,2);
            
            %% Denoise
            CoG = movmean(CoG,3);        
            CoG = CoG - mean(CoG);
            
            %%
            CoGAll{PulseDurationInd,RepetitionInd}(PosInd,:) = CoG;
            
        end
    end
end

%% Plot deltat-position maps & cross-sections
close all;

for PulseDurationInd = 1:PulseListLength

    for RepetitionInd = 1:TotalRepetitions          

        figure(PulseDurationInd)
        subplot(1,TotalRepetitions,RepetitionInd)
        plot(Time,CoGAll{PulseDurationInd,RepetitionInd},'linewidth',2)
        xlabel('Time delay [ps]')
        ylabel('Energy shift [eV]')
        xlim([min(Time) max(Time)])        
        axis square
        ax = gca;
        ax.FontSize = 18;
        ax.LineWidth = 1;
        box on
        
        CoGAllRep(RepetitionInd,:,:) = CoGAll{PulseDurationInd,RepetitionInd};
                
        PV(PulseDurationInd,RepetitionInd) = max(CoGAllRep(RepetitionInd,:,:)) - min(CoGAllRep(RepetitionInd,:,:));

    end
    
    CoGAllRepMean = squeeze(mean(CoGAllRep,1));
    CoGAllRepSTD = squeeze(std(CoGAllRep,0,1));

    figure(PulseDurationInd*1000)
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
    
    colormap jet
    cmap = colormap;
    cmap = cmap(1:round(length(cmap)/PulseListLength):end,:);
%     cmap = flip(cmap,1);

    figure(100)
    hold on
    errorbar(time,CoGAllRepMean + 1.5*(PulseDurationInd-1),CoGAllRepSTD,'linewidth',2,'color',cmap(PulseDurationInd,:))
    hold off
    xlabel('Time delay [ps]')
    ylabel('Energy shift [eV]')
    xlim([min(Time) 2.3])
    ylim([-1 16])
    legend(flip(PulseList))
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    box on
    
end

for PulseDurationInd = 1:PulseListLength
    figure(100)
    hold on
    plot([-0.75 3],1.5*(PulseDurationInd-1)*[1 1],'--k','linewidth',1)
    hold off
end

figure(200)
hold on
errorbar(PulseListNum,mean(PV,2),std(PV,0,2),'xk','linewidth',2)
hold off
xlabel('Laser pulse energy [nJ]')
ylabel('P-V energy shift [eV]')
xlim([0 400])
ylim([0 3.5])
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
box on

%% Plot time scans
close all;

for PulseDurationInd = 1:PulseListLength
    for RepetitionInd = 1:TotalRepetitions          
        
        DataSetAllRep(RepetitionInd,:,:) = squeeze(DataSetCropAll{PulseDurationInd,RepetitionInd}(1,:,:));
        Energy = squeeze(EnergyCropAll{PulseDurationInd,RepetitionInd}(1,:));
        
    end
    
    dataset = squeeze(mean(DataSetAllRep,1));
    dataset = movmean(dataset,3,2);
    dataset = dataset./trapz(Energy,dataset',2)';

    figure(PulseDurationInd)
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

for PulseDurationInd = 1:PulseListLength
    Fs = zeros(TotalRepetitions,length(t));
    for RepetitionInd = 1:TotalRepetitions          
       
        s = CoGAll{PulseDurationInd,RepetitionInd}(1,:);
        s = s./max(s);
        s = [zeros(1,(length(t)+1)/2),s,zeros(1,(length(t)-1)/2-length(s))];

        Fs(RepetitionInd,:) = fftshift(fft(s));
        
    end
    
    Fs = abs(Fs)./max(abs(Fs),[],2);

    FsMean = mean(Fs,1);
    FsSTD = std(Fs,0,1);
    
    figure(300)
    hold on
    errorbar(f*1e-12,FsMean + (PulseDurationInd-1),FsSTD,'linewidth',2,'color',cmap(PulseDurationInd,:))
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

for PulseDurationInd = 1:PulseListLength
    figure(300)
    hold on
    plot([0 5],(PulseDurationInd-1)*[1 1],'--k','linewidth',1)
    hold off
end
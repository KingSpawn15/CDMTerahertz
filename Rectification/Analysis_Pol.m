clc; close all; clear all;

%% Parameters
%\111_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210617061959\ - Used in figures
%\111_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210614231717\ - Focus 0.8ps around t0

% PathMain = 'measurement_data\Polarization\111_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210617061959\';
PathMain = '..\..\..\..\..\..\Desktop\measurement_data\Polarization\111_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210617061959\';

FolderList = dir(PathMain);

% 1-10_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210614115140
% TotalRepetitions = 1;
% Time = -0.05:0.05:0.2;%[ps] %% Time = -0.05 [ps] represents the ZLP calibration point
% Polarization = 0:2:358;

% 111_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210614231717
% TotalRepetitions = 2;
% Time = -0.05:0.05:0.8;%[ps] %% Time = -0.05 [ps] represents the ZLP calibration point
% Polarization = 0:2:358;

% 111_STEM_x1.5k_StaticSpot_10nJ_wTimeScan20210617061959
TotalRepetitions = 1;
Time = -0.05:0.05:3.75;%[ps] %% Time = -0.05 [ps] represents the ZLP calibration point
Polarization = 0:2:178;

LengthTimeVector = length(Time);
LengthPolVector = length(Polarization);

%% Read .dm4 files
Counter = 3;

H = waitbar(0,'Reading .dm4 files...');
for RepetitionInd = 1:TotalRepetitions            
    for PolariztionInd = 1:LengthPolVector
        for TimeInd = 1:LengthTimeVector
        
            FilePath = [PathMain,FolderList(Counter).name];
            saveFlag = false;
            FileName = nan;

            Counter = Counter + 1;

            [DataSet,Energy,Units] = Readdm4_Pol(FilePath);

            %Crop around ZLP area
            if TimeInd == 1
                [~,EnergyZeroInd] = max(DataSet);
                EnergyZero = Energy(EnergyZeroInd);
            end
            
            EnergyHalfWidth = 5;%[eV]
            EnergyMax = EnergyZero + EnergyHalfWidth;
            EnergyMin = EnergyZero - EnergyHalfWidth;

            EnergyCrop = Energy(Energy >= EnergyMin & Energy <= EnergyMax) - EnergyZero;
            DataSetCrop = DataSet(Energy >= EnergyMin & Energy <= EnergyMax);

            EnergyCropAll(RepetitionInd,PolariztionInd,:) = EnergyCrop;    
            DataSetCropAll(RepetitionInd,PolariztionInd,TimeInd,:) = DataSetCrop;
            
            waitbar(Counter/length(FolderList));
            
        end
    end
end

close(H);

%% Analysis
for RepetitionInd = 1:TotalRepetitions            
    for PolariztionInd = 1:LengthPolVector
        for TimeInd = 1:LengthTimeVector
        
            EnergyCrop = squeeze(EnergyCropAll(RepetitionInd,PolariztionInd,:))';        
            DataSetCrop = squeeze(DataSetCropAll(RepetitionInd,PolariztionInd,TimeInd,:))';

            CoG = trapz(EnergyCrop,EnergyCrop.*DataSetCrop,2)./trapz(EnergyCrop,DataSetCrop,2);
            CoGAll(RepetitionInd,PolariztionInd,TimeInd) = CoG;
                      
            %%
            [MaxVal,MaxInd] = max(DataSetCrop);
            
            UpInd = find(DataSetCrop(MaxInd:end)<=0.5*MaxVal,1,'first');
            UpInd = UpInd + MaxInd - 1;
            DownInd = find(DataSetCrop(1:MaxInd)<=0.5*MaxVal,1,'last');
            FWHM = EnergyCrop(UpInd) - EnergyCrop(DownInd);
            FWHMAll(RepetitionInd,PolariztionInd,TimeInd) = FWHM;  
        
        end
    end
end

%% Plot CoG maps
for RepetitionInd = 1:TotalRepetitions            
    for TimeInd = 1:LengthTimeVector
        
        CoGCurrent(:,TimeInd) = CoGAll(RepetitionInd,:,TimeInd);
        
        for p = 1:LengthPolVector/2
            CoGSplit(p,TimeInd) = (CoGCurrent(p,TimeInd) + CoGCurrent(p+LengthPolVector/2,TimeInd))./2;
        end

        pol = Polarization(1:LengthPolVector/2+1);
        
    end

    CoGDenoise = movmean(CoGCurrent,3,1);
    CoGSplitDenoise = movmean(CoGSplit,3,1);

    time = Time(2:end) - 1.05;

    figure(1+10*RepetitionInd)
    imagesc(time,Polarization,CoGDenoise(:,2:end))
    xlabel('Time delay [ps]')
    ylabel('\Lambda/2 angle [deg]')
    xlim([min(time) max(time)])
    ylim([min(Polarization) max(Polarization)])
    colormap jet
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.XTick = -1:0.5:2.5;
    ax.YTick = 0:20:178;
    box on
    caxis([min(CoGDenoise(:)) max(CoGDenoise(:))])
    colorbar

    figure(2+10*RepetitionInd)
    imagesc(time,pol,CoGSplitDenoise(:,2:end))
    xlabel('Time delay [ps]')
    ylabel('\Lambda/2 angle [deg]')
    xlim([min(time) max(time)])
    ylim([min(pol) max(pol)])
    colormap jet
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.XTick = -1:0.5:2.5;
    ax.YTick = 0:20:90;
    box on
    caxis([min(CoGSplitDenoise(:)) max(CoGSplitDenoise(:))])
    colorbar

end
%% Plot FWHM maps
for RepetitionInd = 1:TotalRepetitions            
    for TimeInd = 1:LengthTimeVector

        FWHMCurrent(:,:) = FWHMAll(RepetitionInd,:,:);
        
        for p = 1:LengthPolVector/2
            FWHMSplit(p,TimeInd) = (FWHMCurrent(p,TimeInd) + FWHMCurrent(p+LengthPolVector/2,TimeInd))./2;
        end

        pol = Polarization(1:LengthPolVector/2+1);
        
    end

    FWHMDenoise = movmean(FWHMCurrent,3,1);
    FWHMSplitDenoise = movmean(FWHMSplit,3,1);

    figure(3+10*RepetitionInd)
    imagesc(time,Polarization,FWHMDenoise(:,2:end))
    xlabel('Time delay [ps]')
    ylabel('\Lambda/2 angle [deg]')
    xlim([min(time) max(time)])
    ylim([min(Polarization) max(Polarization)])
    colormap jet
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.XTick = -1:0.5:2.5;
    ax.YTick = 0:20:178;
    box on
    caxis([min(FWHMDenoise(:)) max(FWHMDenoise(:))])
    colorbar

    figure(4+10*RepetitionInd)
    imagesc(time,pol,FWHMSplitDenoise(:,2:end))
    xlabel('Time delay [ps]')
    ylabel('\Lambda/2 angle [deg]')
    xlim([min(time) max(time)])
    ylim([min(pol) max(pol)])
    colormap jet
    axis square
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 1;
    ax.XTick = -1:0.5:2.5;
    ax.YTick = 0:20:90;
    box on
    caxis([min(FWHMSplitDenoise(:)) max(FWHMSplitDenoise(:))])
    colorbar

end

%% Plot individual time scans
%% 10 deg
RepetitionInd = 1;
PolInd = find(Polarization==10);

EnergyCrop = squeeze(EnergyCropAll(RepetitionInd,PolInd,:))';        
DataSetCrop = squeeze(DataSetCropAll(RepetitionInd,PolInd,:,:))';

dataset = movmean(DataSetCrop(:,2:end)',3,1);
dataset = dataset./trapz(EnergyCrop,dataset',1)';

figure(PolInd)
imagesc(EnergyCrop,time,dataset)
xlabel('Energy shift [eV]')
ylabel('Time delay [ps]')
xlim([min(EnergyCrop) max(EnergyCrop)])
ylim([min(time) max(time)])
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:2.5;
box on
caxis([0 max(dataset(:))])
colorbar
RepetitionInd = 1;

tfInd = find(Time>=2,1,'first');
DataSetAll{1} = dataset(1:tfInd-1,:);
PolIndList(1) = PolInd;

%% 56 deg
PolInd = find(Polarization==56);

EnergyCrop = squeeze(EnergyCropAll(RepetitionInd,PolInd,:))';        
DataSetCrop = squeeze(DataSetCropAll(RepetitionInd,PolInd,:,:))';

dataset = movmean(DataSetCrop(:,2:end)',3,1);
dataset = dataset./trapz(EnergyCrop,dataset',1)';

figure(PolInd)
imagesc(EnergyCrop,time,dataset)
xlabel('Energy shift [eV]')
ylabel('Time delay [ps]')
xlim([min(EnergyCrop) max(EnergyCrop)])
ylim([min(time) max(time)])
colormap jet
axis square
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;
ax.YTick = -1:0.5:2.5;
box on
caxis([0 max(dataset(:))])
colorbar

DataSetAll{2} = dataset(1:tfInd-1,:);
PolIndList(2) = PolInd;
TimeDelay_Experiment = Time(2:tfInd);
Energy_Experiment = EnergyCrop;

%% Plot individual time scans
%All
close all;
RepetitionInd = 1;

while true
    for PolInd = 1:1:LengthPolVector

        EnergyCrop = squeeze(EnergyCropAll(RepetitionInd,PolInd,:))';        
        DataSetCrop = squeeze(DataSetCropAll(RepetitionInd,PolInd,:,:))';

%         dataset = movmean(DataSetCrop(:,2:end)',3,1);
        dataset = DataSetCrop(:,2:end);
        dataset = dataset./trapz(EnergyCrop,dataset,1);

        figure(7)
        imagesc(EnergyCrop,time,dataset')
        xlabel('Energy shift [eV]')
        ylabel('Time delay [ps]')
        xlim([min(EnergyCrop) max(EnergyCrop)])
        ylim([min(time) 1.5])
        colormap jet
        axis square
        ax = gca;
        ax.FontSize = 18;
        ax.LineWidth = 1;
        ax.YTick = -1:0.5:2.5;
        box on
        caxis([0 max(dataset(:))])
        colorbar
        annotation('textbox',[0.2, 0.2, 0.1, 0.1],'String',['\theta = ',num2str(Polarization(PolInd)*2),'^o'],...
        'color',[1 1 1],'LineStyle','none','FontSize',18)
        drawnow
        pause(0.3)
        clf
        
    end
end
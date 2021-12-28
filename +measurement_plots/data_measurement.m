function [psi_measurement, e_w_measurement,...
    t_w_measurement] = data_measurement(measurement_polarization)

load('saved_matrices/rectification.mat');

% find index closest to polarization angle
[~, ind] = min( abs( Polarization - measurement_polarization ) );

RepetitionInd = 1;


for PolInd = ind
    
    e_w_measurement = squeeze(EnergyCropAll(RepetitionInd,PolInd,:))';
    DataSetCrop = squeeze(DataSetCropAll(RepetitionInd,PolInd,:,:))';
    
    %         dataset = movmean(DataSetCrop(:,2:end)',3,1);
    dataset = DataSetCrop(:,2:end);
    dataset = dataset./trapz(e_w_measurement,dataset,1);
    
    t_w_measurement = time;
    psi_measurement = dataset.';
    
    
end


function [R1array, yaxis] = forWaterCorrectionCoeffs(patient_numbers)

R1array = zeros(size(patient_numbers,2),1);
yaxis = zeros(size(patient_numbers,2),1);
for patient = 1:length(patient_numbers)
    
    patient_number = patient_numbers(patient)
    if patient_number < 10
        patient_number = ['P00' num2str(patient_number)];
    elseif patient_number<100
        patient_number = ['P0' num2str(patient_number)];
    else
        patient_number = ['P' num2str(patient_number)];
    end

path = pwd;
type = 'GE';
load([path '\Result_MSwcf2\' patient_number type '_M.mat'])

T1MapPre = images{10};
T1MapPost = images{6};
CBVmapPost = images{6};
CBVmapPre = images{10};
%CBVmap = squeeze(CBVmap(:,:,16));

%
R1 = 1000*(1./T1MapPost - 1./T1MapPre);
%R1 = abs(log(1./T1MapPost - 1./T1MapPre));

roiSS = masks_ROIs.SS;
roiWM = masks_ROIs.WM_DSC;
R1array(patient) = mean(mean(R1(roiSS > 0)));
yaxis(patient) = 100*mean(mean((1./CBVmapPost(roiWM > 0) - 1./CBVmapPre(roiWM > 0)))/(R1array(patient)/1000));


end


end


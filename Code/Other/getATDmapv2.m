function [ATDmap,IBATmap,BATmap,BRTmap,sigmap] = getATDmapv2(AIF,path_DSC,Nmeas,Nslices,CmtWM)
%getATDmap Returns Arterial Tissue Delay map and voxel BAT map
%   all values are in data point number (not in sec or msec)
%   *experimental (try to remove section roi dependence)
% Author: Yong Ik Jeong
% Date: 2018-02-06

%% Get brain section ROIs
%roidata = load('..\SAVES\roiset_20180125v2.mat');
global glblTargetPath;
global injectionNum;
global seqType;
tmpdir = dir([glblTargetPath '\ROIs\*ROI_20180125v2.mat']);
if ~isempty(tmpdir)
    roidata = load([glblTargetPath '\ROIs\' tmpdir(1).name]);
    roistack = roidata.roi_stack;
    
    if size(roistack) ~= Nslices
        error('roi mask and dsc Nslices do not match');
    end
    ATDmap = zeros(size(roistack));
    BATmap = zeros(size(roistack));
    IBATmap = zeros(size(roistack));
    BRTmap = zeros(size(roistack));
    %sigmap = cell([length(unique(roistack))-1 2]);
    sigmap = {};
    count = 0;
    
    isvein = 0;
    if exist([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat'],'file')
        % Get vein mask and vein signal
        veindata = load([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat']);
        veinmask = veindata.veinmask;
        %path_LLEPI = [glblTargetPath '\P' sprintf('%03d',injectionNum) '\IR_LL_EPI_PRE'];
        %path_DSC = [glblTargetPath '\P' sprintf('%03d',injectionNum) '\ep2d_perf'];
        if size(veinmask,3) ~= Nslices
            error('vein mask and dsc Nslices do not match');
        end
        vein_signal = zeros(1,Nmeas);
        vein_signal2 = zeros(1,Nmeas);
        for zz = 1:size(veinmask,3)
            dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(zz-1)*Nmeas),(zz*Nmeas),1)));
            for ii = 1:size(veinmask,1)
                for jj = 1:size(veinmask,2)
                    if veinmask(ii,jj,zz)
                        vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),5,'sgolay',3),[1 Nmeas]);
                        vein_signal2 = vein_signal2 + reshape(smooth(dscstack(ii,jj,:),5,'sgolay',1),[1 Nmeas]);
                        %vein_signal = vein_signal + reshape(dscstack(ii,jj,:),[1 Nmeas]);
                        %vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),10),[1 Nmeas]);
                    end
                end
            end
        end
        vein_signal = vein_signal / sum(veinmask(:));
        vein_signal2 = vein_signal2 / sum(veinmask(:));
        [vIBAT, vBAT, vBRT] = findBAT(vein_signal);
        isvein = 1;
        figure;plot(vein_signal);hold on;plot(vein_signal2);
        title(['vein: ' num2str(vBAT) '-' num2str(vBRT)]);   
        
        %[VEIN.Ct,vBRT] = enhanceCorr(AIF,vein_signal,echoTime);
    end
    
    for zz = 1:size(roistack,3)
        dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(zz-1)*Nmeas),(zz*Nmeas),1)));
        tmproimask = roistack(:,:,zz);
        tmpbatmap = zeros(size(tmproimask));
        tmpibatmap = zeros(size(tmproimask));
        tmpbrtmap = zeros(size(tmproimask));
        tmpatdmap = zeros(size(tmproimask));
        mask = automaskns(imfilter(dscstack,fspecial('gaussian',[3 3],2)));
        sdscstack = fSpatialFilter(dscstack,[2 3],mask);
        
        for ii = 1:size(roistack,1)
            for jj = 1:size(roistack,2)
                if mask(ii,jj)
                    
                    signal = reshape(smooth(sdscstack(ii,jj,:),5,'sgolay',3),[1 Nmeas]);
                    signal2 = reshape(smooth(sdscstack(ii,jj,:),5,'sgolay',1),[1 Nmeas]);
                                        
                    [ibat,bat,brt] = findBAT(signal);
                    %[~,brt] = enhanceCorr(AIF,signal,echoTime);
                    if brt == 0 || ibat == 0 || bat == 0
                        1;
                    end
%                     if ismember(nonzero(ind),[7, 11, 23, 28])
%                         1;
%                     end
                    % 20180305 YIJ: if bat and brt are greater than vein, set
                    % to average wm bat and brt...
                    if isvein
                        if bat > vBAT + 5 % 20180307 YIJ: +5 experimental
                            if roistack(ii,jj,zz) == 7
                                1;
                            end
%                             figure;plot(signal);hold on;plot(signal2);
%                             %title(['ROI: ' num2str(nonzero(ind)) ' bat ' num2str(bat) '-' num2str(vBAT)]);
%                             title(['bat ' num2str(bat) '-' num2str(vBAT)]);
                            
                            count = count + 1;
                            sigmap{count,1} = count;
                            sigmap{count,2} = roistack(ii,jj,zz);
                            sigmap{count,3} = signal;
                            sigmap{count,4} = [ibat bat brt vIBAT vBAT vBRT];
                            %bat = CmtWM(2);
                            bat = vBAT;
                        end
                        if brt > vBRT + 5
                            if roistack(ii,jj,zz) == 7
                                1;
                            end
%                             figure;plot(signal);hold on;plot(signal2);
%                             %title(['ROI: ' num2str(nonzero(ind)) ' brt ' num2str(brt) '-' num2str(vBRT)]);
%                             title(['brt ' num2str(brt) '-' num2str(vBRT)]);
%                             
                            count = count + 1;
                            sigmap{count,1} = count;
                            sigmap{count,2} = roistack(ii,jj,zz);
                            sigmap{count,3} = signal;
                            sigmap{count,4} = [ibat bat brt vIBAT vBAT vBRT];
                            %brt = CmtWM(3);
                            brt = vBRT;
                        end
                    end
                    tmpbatmap(ii,jj) = bat;
                    tmpibatmap(ii,jj) = ibat;
                    tmpbrtmap(ii,jj) = brt;
                    tmpatdmap(ii,jj) = bat - AIF.BATP;
%                     sigmap{count,1} = count;
%                     %sigmap{count,2} = nonzero(ind);
%                     sigmap{count,2} = roistack(ii,jj,zz);
%                     sigmap{count,3} = signal;
                end
            end
        end
        ATDmap(:,:,zz) = tmpatdmap;
        BATmap(:,:,zz) = tmpbatmap;
        IBATmap(:,:,zz) = tmpibatmap;
        BRTmap(:,:,zz) = tmpbrtmap;
    end
else
    ATDmap = [];
    BATmap = [];
    IBATmap = [];
    BRTmap = [];
    sigmap = [];
end

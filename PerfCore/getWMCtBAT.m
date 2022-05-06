function [ wm_signal, CmtWM, nslice ] = getWMCtBAT( path_DSC,N_slices,N_meas )
%getWMCtBAT Get bolus arrival times and Ct curve of white matter
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-08-14

global glblTargetPath;
global injectionNum;
global seqType;
if exist([glblTargetPath '\WM_Mask_P' sprintf('%03d',injectionNum) seqType '.mat'],'file')
    wmdata = load([glblTargetPath '\WM_Mask_P' sprintf('%03d',injectionNum) seqType '.mat']);
    WM_SS = wmdata.WM_SS;
    aifdata = load([glblTargetPath '\AIFdata\P' sprintf('%03d',injectionNum) seqType(1:2) '_AIF_A.mat']);
    %tmpind = strfind(path_DSC,'\');
    path_LLEPI = [glblTargetPath '\P' sprintf('%03d',injectionNum) '\IR_LL_EPI_PRE'];
    if isfield(wmdata,'nslice')
        nslice = wmdata.nslice;
    else
        [mask, nslice, match_index] = fAutoMaskWM_DSC_Philips(path_LLEPI,path_DSC,WM_SS,N_slices,N_meas);
    end
    dscslice = aifdata.ROIs.dsc_stack{nslice};
    wm_signal = zeros(1,size(dscslice,3));
    for ii = 1:size(dscslice,1)
        for jj = 1:size(dscslice,2)
            if WM_SS(ii,jj)
                %wm_signal = wm_signal + freqfilter(reshape(smooth(dscslice(ii,jj,:),5,'sgolay',3),1,size(dscslice,3)),'low',10);
                wm_signal = wm_signal + reshape(smooth(dscslice(ii,jj,:),.05,'sgolay',3),[1 size(dscslice,3)]);
            end
        end
    end
    wm_signal = wm_signal / sum(WM_SS(:));
    %tmp_wm_signal = freqfilter(wm_signal,'low',10);
    [ibat, bat, brt] = findBAT(wm_signal);
    %YIJ 20170829 try use peak location instead of brt
    [~,nmin] = min(wm_signal(:));
    CmtWM = [ibat bat brt nmin];
    aifbat = aifdata.cutoffs_ROIs.AIF;
    aifwindow = aifbat(3)-aifbat(2);
    wmwindow = CmtWM(3)-CmtWM(2);
    if aifwindow > wmwindow
        CmtWM(3) = CmtWM(3) + (aifwindow - wmwindow);
    end
else
    wm_signal = [];
    CmtWM = [];
    nslice = [];
end

    
end


function [noDD,DD,CmtWM,etc,sigmap,VEINslice,veinmask,positionVein,n_Veinslice] = Calc_nSVD_Philips_DD_all(path_DSC,N_slices,signal_AIF,cutoffs_AIF,N_meas,AIFslice,n_slice_AIF)

%%=======================================================================%%
%%          rCBV,rCBF,MTT,1/SNR_c calculation function         
%%=======================================================================%%
% input :   1) path_DSC (string)  : path of DSC images
%           2) N_slices (scalar)  : total slice number in DSC images
%           3) signalAIF (1xN vector) : Signal intensity for AIF
%           4) cutoffs_AIF ([i,j]) : BATP and RTP in AIF signal
%%=======================================================================%%
% output:   1) CBF_SVD  (3-D martrix) : rCBF in SVD
%           2) CMMT_SVD (3-D martrix) : CMTT from Residue fuction in SVD
%           3) CBV_CVP  (3-D martrix) : = CBF_SVD * CMMT_SVD 
%           4) CBV_DSC  (3-D martrix) : rCBV from the ratio of integral of
%           5) CMTT_CVP (3-D martrix) : = CBV_SDSC / CBF_SVD 
%           6) reciSNR_c(3-D martrix) : 3-D SNR in Conct.
%                                                               02-22-2004 

% read header of first image
% header = dicominfo([path_DSC '\1.dcm']);
header = dicominfo([path_DSC '\1.dcm'],'Dictionary','dicom-dict.txt'); %Upon the upgrade of the dicom dictionary
Dt = (header.RepetitionTime)/1000; % from ms to sec
echoTime = header.EchoTime;
Wij_SVD = 0.2;
%disp('1');
% resize image
tempimage = dicomread([path_DSC '\1.dcm']);
%resize in the 3 cases where images are 112 %%%%Grady
% if (size(tempimage,1)<=112) || (size(tempimage,2)<=112)
%     tempimage = imresize(tempimage,128/112);
% else
%     tempimage = tempimage;
% end
    
global sampdelayfix;
global DDfix;
global MANUAL;
if sampdelayfix
    fprintf('calc svd sampdelayfix\n');
    AIF_struct = Calc_AIFstruct(path_DSC,AIFslice,n_slice_AIF,N_slices);
    AIF = AIFslice;
else
% calculate AIF struct
AIF = Sig2Conct_AIF(signal_AIF,cutoffs_AIF, echoTime); % to Wanyong's code
end

if DDfix
    fprintf('DD fix on\n');
else
    fprintf('DD fix off\n');
end

% definition for defragment

noDD.CBV_DSC  = zeros(size(tempimage,1),size(tempimage,2),N_slices);
noDD.CBF_SVD = zeros(size(tempimage,1),size(tempimage,2),N_slices); 
noDD.CMTT_CVP = zeros(size(tempimage,1),size(tempimage,2),N_slices);
noDD.Tmax_map = zeros(size(tempimage,1),size(tempimage,2),N_slices); 
noDD.dBATmap = zeros(size(tempimage,1),size(tempimage,2),N_slices); 
noDD.TTP_map = zeros(size(tempimage,1),size(tempimage,2),N_slices);
noDD.isbolus = zeros(size(tempimage,1),size(tempimage,2),N_slices);
DD = [];
if DDfix
DD.CBV_DSC  = zeros(size(tempimage,1),size(tempimage,2),N_slices);
DD.CBF_SVD = zeros(size(tempimage,1),size(tempimage,2),N_slices); 
DD.CMTT_CVP = zeros(size(tempimage,1),size(tempimage,2),N_slices);
end


ATDmap = [];
BATmap = [];
IBATmap = [];
BRTmap = [];
sigmap = [];
VEIN = [];
VEINslice = [];
veinmask = [];
positionVein = [];
n_Veinslice = [];

[wm_signal, CmtWM, t1slc] = getWMCtBAT(path_DSC,N_slices,N_meas);
if isempty(CmtWM)
    if MANUAL
        error('Cant find WM mask for finding WM Ct and BAT');
    else
        warning('Cant find WM mask for finding WM Ct and BAT');
    end
end

dscstacks = [];
for k = 1:N_slices
    tempimages = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(k-1)*N_meas),(k*N_meas),1))); %80 measurements
    images = tempimages;
    
    %YIJ 20170615 smooth mask
    mask = automaskns(imfilter(images,fspecial('gaussian',[3 3],2)));
    
    % 20210427 YIJ spatial smooth is now done in smoothDSC.m
    % spatial filtering
    dscstacks{k} = images;%fSpatialFilter(images,[2 3],mask);
    automasks{k} = mask;
end

% 20210326 YIJ: spatially smooth DSC (xy and z)
dscstacks = smoothDSC(dscstacks);

if DDfix
    [VEINslice,veinmask,positionVein,n_Veinslice] = getVEIN(Dt,echoTime,path_DSC,N_meas,N_slices,t1slc,n_slice_AIF,AIFslice);
    VEIN_struct = Calc_AIFstruct(path_DSC,VEINslice,n_Veinslice,N_slices);
    VEIN = VEINslice;
    if isempty(t1slc)
        t1slc = n_Veinslice;
    end
end

    [ATDmap,IBATmap,BATmap,BRTmap,sigmap] = getATDmap_test(AIF,dscstacks,automasks,path_DSC,N_meas,N_slices,CmtWM,t1slc,Dt,VEIN);
    if isempty(ATDmap)
        fprintf('**Tried to do DDfix but no ATDmap could be found\n');
    end
    etc.ATDmap = ATDmap;
    etc.IBATmap = IBATmap;
    etc.BATmap = BATmap;
    etc.BRTmap = BRTmap;

% calculate & show remaining time 
%handle_progress = ProgressBar('Calculating normal SVD...');

%tmpaif = load('..\saves\average global aif_20171208.mat');

%Jessy (4/29/2009):
%N_meas = 200;
% just analysis of a slice on WM_DSC ROI
% if DDfix
%     [VEIN] = getVEIN(Dt,echoTime,path_DSC,N_meas,N_slices,t1slc);
% end
if DDfix
            % get delay and dispersion fit using resampled AIF
            %[DDfit,VEIN] = DDcorr_test(AIF,Dt,echoTime,Wij_SVD,path_DSC,N_meas,N_slices,t1slc);
            VEIN.Ct = VEIN_struct{t1slc};
            AIF.Ct = AIF_struct{t1slc};
            [DDfit] = DDcorr_test(AIF,Dt,Wij_SVD,VEIN,path_DSC);
            etc.DDfit = DDfit;
            if isempty(DDfit)
                fprintf('**Tried to do DDfix but no fit could be found\n');
            end
else
    DDfit = [];
end
        
for k = 1:N_slices %.................................
    
    %DDfit = [];
    if sampdelayfix
        AIF.Ct = AIF_struct{k};
        if ~isempty(ATDmap)
            kATDmap = ATDmap(:,:,k);
            kBATmap = BATmap(:,:,k);
            kIBATmap = IBATmap(:,:,k);
            kBRTmap = BRTmap(:,:,k);
        end
        
%         if DDfix
%             % get delay and dispersion fit using resampled AIF
%             %[DDfit,VEIN] = DDcorr_test(AIF,Dt,echoTime,Wij_SVD,path_DSC,N_meas,N_slices,t1slc);
%             VEIN.Ct = VEIN_struct{k};
%             [DDfit] = DDcorr_test(AIF,Dt,Wij_SVD,VEIN);
%             if isempty(DDfit)
%                 fprintf('**Tried to do DDfix but no fit could be found\n');
%             end
%         end
        
    end
    %AIF = tmpaif.AIF;
    
    if k== 15 %20190718 YIJ: test for case 23
        1;
    end
    
%     % read images fast
%      %tempimages = double(Cells2Matrix(IA_ReadImages(path_DSC,k,0,N_slices)));
%     tempimages = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(k-1)*N_meas),(k*N_meas),1))); %80 measurements
%     
% %     if (size(tempimages,1)<=112) || (size(tempimages,2)<=112)
% %         images = imresize(tempimages,128/112);
% %     else
%         images = tempimages;
% %     end
%     
%     % make a non-scalp mask
%     %mask = automaskns(images); % to Wanyong's code
%     
%     %YIJ 20170615 smooth mask
%     mask = automaskns(imfilter(images,fspecial('gaussian',[3 3],2)));
%     
%     % spatial filtering
%     S_images = fSpatialFilter(images,[2 3],mask);
%     % 20190725 YIJ: Try filter twice (might be better for Tmax)
%     %S_images = fSpatialFilter(S_images,[2 3],mask);
%     %S_images = images;

    S_images = dscstacks{k};
    mask = automasks{k};

    % calculate MTT_SVD
    if k == 2
        1;
    end

    if DDfix && (~isempty(ATDmap) && ~isempty(DDfit))
        %debugDDaif(S_images,AIF,Dt,Wij_SVD,echoTime,DDfit,kATDmap,kIBATmap,kBATmap,kBRTmap,k,N_slices);
    end
    
    if ~isempty(ATDmap)
        if DDfix && isempty(DDfit)
            %YIJtmp error('DDfix is on, but DDfit is empty');
        end
        [tmpnoDD,tmpDD] = fCBV_CBF_MTT_nSVD_DD_all(S_images,mask,AIF,Dt,Wij_SVD,echoTime,CmtWM,DDfit,kATDmap,kIBATmap,kBATmap,kBRTmap);
    else
        error('Did not run fCBV_CBF_MTT_nSVD_DD_all');
    end
    
    
    noDD.CBV_DSC(:,:,k)   = tmpnoDD.rCBV_DSC;
    noDD.CBF_SVD(:,:,k)   = tmpnoDD.rCBF_SVD;
    noDD.CMTT_CVP(:,:,k)  = tmpnoDD.CMTT_CVP;
    noDD.Tmax_map(:,:,k)   = tmpnoDD.Tmax_map;
    noDD.dBATmap(:,:,k)   = tmpnoDD.dBATmap;
    noDD.TTP_map(:,:,k) = tmpnoDD.TTP_map;
    noDD.isbolus(:,:,k) = tmpnoDD.isbolus;
    if DDfix
    DD.CBV_DSC(:,:,k)   = tmpDD.rCBV_DSC;
    DD.CBF_SVD(:,:,k)   = tmpDD.rCBF_SVD;
    DD.CMTT_CVP(:,:,k)  = tmpDD.CMTT_CVP;
    end

    
    %ProgressBar(handle_progress,k/N_slices)
end

%close(handle_progress)

% clear variables
clear header Dt Wij_SVD AIF mask images Wf S_images temp_CBV_DSC temp_CBF_SVD temp_CMTT_CVP tempTmaxMap temp_dBATmap
 
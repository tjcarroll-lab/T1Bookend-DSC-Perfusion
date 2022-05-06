function [Signal_AIF,cutoffs_AIF,position,Conct_AIF,Images,n_AIFslice,AIF,positionAIF,Mask_AIF_slice] = autoaif_wy_Philips_v2(path_DSC,N_Slices,N_meas,tmpdata)

% [Signal_AIF, cutoffs_AIF, position, Conct_AIF] = autoaif_wy(path_DSC)
%
% Automated optimal AIF algorithm 
% Author: Wanyong
%
% Description: aterial input function is chosen from automated algorithm 
% based on the following categories satisfied 
% A. exclude the voxels to satify
% 1. remove saturation 
% 2_1. remove late BAT based on Maximum concentration
% 2_2. remove late BAT based on BAT
% 3. remove late recirculation
% 4. remove smaller integral of Conct. up to recir.
% 5. remove lower highest peak.
% 6. remove wider bolus.1. the normalized 4-point integral over time (starting at the voxel 
% B. choose the more than 5 voxels to satisfy the following
% 1.
% 
% Status: stable

% Versions
%   06-06-29 (WS): initial version
%   07-01-10 (WS): Modified
%--------------------------------------------------------------------------

global glblTargetPath;
global injectionNum;
global seqType;
global MANUAL;
global caseid;
global currtime;

% caseid = split(path_DSC,{'\','/'});
%     caseid = caseid{end-3};
%     currtime = datestr(now,'yyyymmddHHMMSS');
    
if exist([glblTargetPath '\AIF_Mask_P' sprintf('%03d',injectionNum) seqType '.mat'],'file')
  
    aifdata = load([glblTargetPath '\AIF_Mask_P' sprintf('%03d',injectionNum) seqType '.mat']);
    aifmask = aifdata.aifmask;
    
    if size(aifmask,3) ~= 1 && size(aifmask,3) ~= N_Slices
        error('vein mask and dsc Nslices do not match');
    end
    if size(aifmask,3) ~= 1
        [tmpi,tmpj,tmpz] = ind2sub(size(aifmask),find(aifmask,1));
        %Mask_AIF_slice = aifmask(:,:,tmpz);
        Mask_AIF_slice = aifmask;
        n_AIFslice = tmpz;
    else
        error('aif mask is single slice');
    end
    
    Signal_AIF = zeros(1,N_meas);
    Signal_AIF2 = zeros(1,N_meas);
    Conct_AIF = zeros(1,N_meas);
    
    for s = 1:N_Slices
%     tempimage = double(Cells2Matrix(IA_ReadImages(path_DSC,s,0,N_Slices)));
        tempimage = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(s-1)*N_meas),(s*N_meas),1))); 
        Images{s} = tempimage;
    end
    
    % Spatially smooth DSC (xy and z)
    Images = smoothDSC(Images);
    
    header = dicominfo([path_DSC '\1.dcm'],'Dictionary','dicom-dict.txt');
    TR = header.RepetitionTime; %ms    
    TE = header.EchoTime;  %ms

    dscstack = Images{n_AIFslice};
    for ii = 1:size(Mask_AIF_slice,1)
        for jj = 1:size(Mask_AIF_slice,2)
            if Mask_AIF_slice(ii,jj,n_AIFslice)
                                
                % 20210422 YIJ: remove smoothing because of arrival time shift
                tmpsignal = reshape(dscstack(ii,jj,:),[1 N_meas]);
                Signal_AIF = Signal_AIF + tmpsignal;
                % 20210329 YIJ: do spline smoothing
%                 tmpfit = fit([1:N_meas]',squeeze(dscstack(ii,jj,:)),'smoothingspline','smoothingparam',.8);
%                 vein_signal = vein_signal + reshape(tmpfit(1:N_meas),[1 N_meas]);
%                 vein_signal(vein_signal < 0) = 0;
                Signal_AIF2 = Signal_AIF;
                
                [ibat,bat,brt] = findBAT_test(tmpsignal);
                tmpconct = -1/TE*log(tmpsignal./mean(tmpsignal(ibat:bat))); % temp = conct.
                tmpconct(find(~isfinite(tmpconct)))=0;
                Conct_AIF = Conct_AIF + tmpconct;
            end
        end
    end
    
    Signal_AIF = Signal_AIF/sum(Mask_AIF_slice(:));
    Signal_AIF2 = Signal_AIF2/sum(Mask_AIF_slice(:));
    Conct_AIF(Conct_AIF < 0) = 0;
    Conct_AIF = Conct_AIF/sum(Mask_AIF_slice(:));
    
    [ibat bat brt] = findBAT_test(Signal_AIF2);
    
    % YIJ 20180910
    if bat < 1 || brt >= 199
        figure;plot(Signal_AIF2);title(['bad aif on slc ' num2str(n_AIFslice)]);
        error('bad aif');
    end
    cutoffs_AIF = [ibat bat brt];
    
    figure;plot(Signal_AIF);hold on;plot(Signal_AIF2);
    title(['aif: ' num2str(bat) '-' num2str(brt) ' TR: ' num2str(TR/1000)]);
    %exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_aif(manu).jpg'], 'resolution', 300);
    
    count = 1;
    clear positionAIF
    for i = 1:size(Mask_AIF_slice,1)
        for j = 1:size(Mask_AIF_slice,2)
            if Mask_AIF_slice(i,j,n_AIFslice)
                positionAIF{count} = [i,j,n_AIFslice];
                count = count + 1;
            end
        end
    end
    position = positionAIF;
    
    Dt = TR*1e-3;
    AIF = Conct_AIF(bat:brt);
    time = 0:Dt:(size(AIF,2)-1)*Dt;
    
    time2 = time;
    AIF2 = AIF;
    
    BestFit = gamma_variate_IDC_gridsearch(time2,AIF2);
    SumSquares = BestFit(1);
    A = BestFit(2);
    alpha = BestFit(3);
    beta = BestFit(4);
    [A,alpha,beta,RESNORM,RESIDUAL,EXITFLAG] = gamma_variate_IDC_LeastSquaresFit(time2,AIF2,A,alpha,beta);
    FitParam.A = A;
    FitParam.alpha = alpha;
    FitParam.beta = beta;
    time_smooth = 0:Dt/N_Slices:(brt-bat)*Dt;
    SmoothAIFfit = gamma_variate_model_IDC(time_smooth,A,alpha,beta);
    figure();
    plot(time2,AIF2,'r*',time_smooth,SmoothAIFfit,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('AIF - Smooth Fit');
    legend('Measured AIF','Gamma variate fit');
    %exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_aifFit(manu).jpg'], 'resolution', 300);
    
    %yv(dscstack(:,:,n_AIFslice), [], 'ol', Mask_AIF_slice(:,:,n_AIFslice))    
    figure;imshow(dscstack(:,:,n_AIFslice),[],'initialmagnification','fit');
    roir = ones(size(dscstack,[1 2]));
    roib = zeros(size(dscstack,[1 2]));
    roig = zeros(size(dscstack,[1 2]));
    roic = cat(3, roir, roig, roib);
    hold on; imh = imshow(roic);
    set(imh, 'alphadata', 0.5*Mask_AIF_slice(:,:,n_AIFslice));
    
%exportgraphics(gca, ['..\images\autovein debug\' currtime '_' caseid '_aifloc(manu).jpg'], 'resolution', 300);

    area = sum(SmoothAIFfit);
    AIF = struct('Ct',SmoothAIFfit,'area',area,'FitParam',FitParam,'RESNORM',RESNORM,'ITP',ibat,'BATP',bat,'RTP',brt,'signal',Signal_AIF,'conc',Conct_AIF);

else

n_AIFslice = 0;
AIF = [];
positionAIF = [];
flg_plot = 0;

% % read basic info.
% [DSCimage, header, N_Slices] = ReadDICOMSeries(path_DSC,1);
header = dicominfo([path_DSC '\1.dcm'],'Dictionary','dicom-dict.txt'); %Upon the upgrade of the dicom dictionary


%basic info for anonimized images
%N_Slices = 5;
% N_Slices = 15;
%N_meas = 200;
for  s = 1:N_Slices
%     tempimage = double(Cells2Matrix(IA_ReadImages(path_DSC,s,0,N_Slices)));
    tempimage = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(s-1)*N_meas),(s*N_meas),1))); 
    %DO NOT RESIZE YET because resolution = 96 x 96
%     if (size(tempimage,1)>=256) || (size(tempimage,2)>=256)
%         Images{s} = imresize(tempimage,0.5);
%     else
        Images{s} = tempimage;
%     end
end

% 20210326 YIJ: spatially smooth DSC (xy and z)
Images = smoothDSC(Images);

%Get TR and TE
% TR = header(1).RepetitionTime; %ms    
% TE = header(1).EchoTime;       %ms
TR = header.RepetitionTime; %ms    
TE = header.EchoTime;       %ms
timepoints = size(Images{1},3);


% find BAT
NmaskNS = 0;
%f1 = figure;
for s = 1:N_Slices
    % read images
    images = Images{s};
    %Create mask 
    %YIJ 20170614 smooth mask
    %mask_s = automaskns(images,0.1);
    mask_s = automaskns(imfilter(images(:,:,1),fspecial('gaussian',[3 3],2)));
    %mask_s = images(:,:,1) > 0;
    
    % YIJ 20170509:
    % Make mask without scalp manually
    %figure(f1); imshow(images(:,:,1), []); title('Make mask of brain for aif'); mask_s = roipoly; mask_all{s} = mask_s;
    
    sum_signal = zeros(1,size(images,3));
    for i = 1:size(images,1)
    for j = 1:size(images,2)
    if mask_s(i,j)
        %sum_signal = sum_signal+ reshape(images(i,j,:),[1 size(images,3)]);
        %tmpary = reshape(images(i,j,:),[1 size(images,3)]);
        %YIJ 20170614 smooth & take out possible leakage/T1 enhance
        
        %tmpary = reshape(smooth(images(i,j,:),max(sqrt(1/(TR/1000))*5,5),'sgolay',3),[1 size(images,3)]);
        
        %tmpary = freqfilter(tmpary,'low',10);
%         if mean(tmpary(1:timepoints/10)) > 0.95*mean(tmpary(end-timepoints/10:end))
%             sum_signal = sum_signal+ tmpary;
%         else
%             mask_s(i,j) = 0;
%         end

        % 20210422 YIJ: remove smoothing because of arrival time shift
        tmpary = reshape(images(i,j,:),[1 timepoints]);
        % 20210329 YIJ: do spline smoothing
%         tmpfit = fit([1:timepoints]',squeeze(images(i,j,:)),'smoothingspline','smoothingparam',.8);
%         tmpary = reshape(tmpfit(1:timepoints),[1 timepoints]);
%         tmpary(tmpary < 0) = 0;
            
        sum_signal = sum_signal+ tmpary;
    end
    end
    end
%     size(sum_signal)
    Sum_Signal(s,:) =  sum_signal;
    NmaskNS = NmaskNS + sum(sum(mask_s));
end
%close(f1);
Mean_Signal = sum(Sum_Signal)/NmaskNS;

%Time axis
time = 0:TR*1e-3:(size(images,3)-1)*TR*1e-3;

% pre definition
TimePoint           = size(images,3);
MaxCtmap            = zeros(size(images,1),size(images,2),N_Slices);
% IntegralCt2BRTmap   = zeros(size(images,1),size(images,2),N_Slices);
AreaoverHightmap    = zeros(size(images,1),size(images,2),N_Slices);
% StdsCt2Maxmap       = zeros(size(images,1),size(images,2),N_Slices);
Mask_AIF            = zeros(size(images,1),size(images,2),N_Slices);

cond_map = zeros(size(images,1),size(images,2),N_Slices,10);
tmpneg_map = zeros(size(images,1),size(images,2),N_Slices);

% find BAT from Mean Signal
[iBAT BAT BRT] = findBAT_test(Mean_Signal);
if(iBAT == 0)
    iBAT = 1;
end

%figure();plot(Mean_Signal);title([num2str(iBAT) ' ' num2str(BAT) ' ' num2str(BRT)]);

[Smin Nmin] = min(Mean_Signal);
[Smax Nmax] = max(Mean_Signal);

% convert Gd_concent
MeanS0 = mean(Mean_Signal(iBAT:BAT));
Mean_Conct = -1/TE*log(Mean_Signal./MeanS0); 
Mean_Conct(find(~isfinite(Mean_Conct)))=0;
%YIJ 20170608
Mean_Conct(Mean_Conct < 0) = 0;

MeanCtmax  = max(Mean_Conct);

warning off

% mean chage in conct. in a whole brain
Mean_Gd = -1/TE*log(Mean_Signal./mean(Mean_Signal(iBAT:BAT)));
Mean_Gd(find(~isfinite(Mean_Gd)))=0;
%      f1 = figure;
%      f2 = figure;
% make Conct. map
for s = 1:N_Slices
    % read images
    images = Images{s};

    %Create new mask 
    %YIJ 20170614 smooth mask
    %mask_AIF = automaskns(images(:,:,1),0.1);
    mask_AIF = automaskns(imfilter(images(:,:,1),fspecial('gaussian',[3 3],2)));
    %mask_AIF = images(:,:,1) > 0;
    
    % YIJ 20170509: use same mask as above
    %mask_AIF = mask_all{s};
    %figure(); imshow(images(:,:,BAT+10), []); title('Make mask of brain'); mask_AIF = roipoly;
    
    Gd_con{s} = zeros(size(images));
    allbat{s} = zeros(size(images));
    filt_sig{s} = zeros(size(images));
    
    %YIJ 20170616 parallel computing
    jsize = size(images,2);
    tmpgd_con = zeros(size(images,1),jsize,size(images,3));
    tmpbat = zeros(size(images,1),jsize,3);
    tmpfilt_sig = zeros(size(images,1),jsize,size(images,3));

%     if ~matlabpool('size')
%         matlabpool
%     end
    % Convert signals in all voxels to [Gd] using mask
    %parfor i = 1:size(images,1)
    for i = 1:size(images,1)
    for j = 1:jsize
        if i == 180 && j == 130 && s == 2
            %disp('debug');
        end
        if mask_AIF(i,j)
            %temp = reshape(images(i,j,:),[1 TimePoint]);
            %YIJ 20170614 smooth
            %temp = reshape(smooth(images(i,j,:),5,'sgolay',3),[1 TimePoint]);
            %temp2 = reshape(smooth(images(i,j,:),max(sqrt(1/(TR/1000))*5,5),'sgolay',3),[1 TimePoint]);
            
            % 20210422 YIJ: remove smoothing because of arrival time shift
            temp = reshape(images(i,j,:),[1 timepoints]);
            % 20210329 YIJ: do spline smoothing
%             tmpfit = fit([1:timepoints]',squeeze(images(i,j,:)),'smoothingspline','smoothingparam',.8);
%             temp = reshape(tmpfit(1:timepoints),[1 timepoints]);
            temp2 = temp;
            
            %temp2 = temp;
            %YIJ 20170626 smoothing can cause negative numbers sometimes
            temp(temp < 0) = 0;
            temp2(temp2 < 0) = 0;
            
            %temp = freqfilter(temp,'low',10);
%             tmpfilt_sig(i,j,:) = temp;
%             if sum(temp < 0)
%                 tmpneg_map(i,j,s) = 1;
%                 temp(temp < 0) = 0;
%             end
%             tmptmp = temp;
%             tmptmp2 = smooth(temp,5,'sgolay',3)';
%             [ibat2 bat2 brt2] = findBAT(tmptmp2);
            
%             %YIJ 20170626 take out possible leakage/T1 enhance
%             if mean(temp2(1:timepoints/10)) > 0.95*mean(temp2(end-timepoints/10:end))
%                 %figure(1);plot(temp);pause;
%             else
%                 mask_AIF(i,j) = 0;
%             end
            
            meanS0 = mean(temp(iBAT:BAT));
            stdS0 = std(temp(iBAT:BAT));
            
%             if meanS0 < 400
%                 mask_AIF(i,j) = 0;
%             end
                

            [smax nmax] = max(temp);   % temp = signal
            [smin nmin] = min(temp);   % temp = signal

            % calcualte iBAT, BAT & BRT
            [ibat bat brt] = findBAT_test(temp2);
            
            if ibat & bat & brt
            
                % Low baseline signal is probably wrong
                if mean(temp(ibat:bat)) < 200
                    mask_AIF(i,j) = 0;
                end
                
                % make Gd_concent maps
                temp = -1/TE*log(temp./mean(temp(ibat:bat))); % temp = conct.
                temp(find(~isfinite(temp)))=0;
                tmpgd_con(i,j,:) = temp;
                %YIJ 20170608
                temp(temp < 0) = 0;
                
                temp2 = -1/TE*log(temp2./mean(temp2(ibat:bat))); % temp = conct.
                temp2(find(~isfinite(temp2)))=0;
                %YIJ 20170608
                temp2(temp2 < 0) = 0;
                
                %temp2 = temp;
                
                
                tmpbat(i,j,:) = [ibat bat brt];
            
                % save data
                Ctmax                       = max(temp);
                MaxCtmap(i,j,s)             = Ctmax;
%                 IntegralCt2BRTmap(i,j,s)    = sum(temp(bat:brt));
                AreaoverHightmap(i,j,s)     = sum(temp(bat:brt))/Ctmax;
%                 StdsCt2Maxmap(i,j,s)        = std(temp(iBAT:BAT-1));
            
                % 1. remove saturation 
                %if smin < smax*0.05+ 2.33*stdS0 % 99% of possibility
                if smin < smax*0.1+ 2.33*stdS0    
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,1) = 1;
%                     [pks00,locs00] = findpeaks(temp,'minpeakdistance',round(timepoints/20),'minpeakheight',0.6*Ctmax);
% %                     if length(pks00) == 1
%                     [pks0,locs0] = findpeaks(temp(max([nmin-3 1]):min([nmin+3 TimePoint])),'minpeakheight',0.6*Ctmax);
%                     %[pks0,locs0] = findpeaks(temp,'minpeakheight',0.6*Ctmax);
%                     if ((nmin > bat) && (nmin < brt)) %&& (sum((locs0 > nmin-3) & (locs0 < nmin+3)) == 1)
%                     if length(pks00) == 1 && length(pks0) == 1
%                         figure(f1);
%                         plot(temp);
%                         figure(f2);
%                         plot(reshape(smooth(images(i,j,:),5,'sgolay',3),[1 TimePoint]));
%                         title([num2str(tmpbat(i,j,:)) ' ' num2str([i j s])]);
%                         pause;
%                     end
%                     end
                    %figure(2);plot(tmptmp);title([num2str(ibat) ' ' num2str(bat) ' ' num2str(brt)]);pause;
                end
                % 2_1. remove late BAT based on Maximum concentration
                if nmin > Nmin 
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,2) = 1;
                end
                % 2_2. remove late BAT based on BAT
                if bat > BAT + 1 % exclude more than 2 secs of BAT
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,3) = 1;
                end
                % 3. remove late recirculation
                if brt > BRT + 2 % exclude more than 3 secs of BRT 
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,4) = 1;
                end
                % 4. remove smaller integral of Conct. up to recir.
                if sum(temp(bat:brt)) < sum(Mean_Conct(BAT:BRT))
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,5) = 1;
                end
                % 5. remove lower highest peak.
                if Ctmax < max(Mean_Conct)
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,6) = 1;
                end
                % 6. remove wider bolus.
                if sum(temp(bat:brt))/Ctmax > sum(Mean_Conct(BAT:BRT))/MeanCtmax
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,7) = 1;
                end
                
                %YIJ 20170818 bolus peak should be between bat and brt
                if nmin < bat || nmin > brt
                    mask_AIF(i,j) = 0;
                    cond_map(i,j,s,8) = 1;
%                 else
%                     try
%                     [pks2,locs] = findpeaks(temp(bat:nmin+1),'minpeakdistance',round(length(bat:nmin+1)/10));
%                     catch
%                         disp('hi');
%                     end
%                     if length(pks2) > 1
%                         figure(f1);plot(reshape(smooth(images(i,j,:),5,'sgolay',3),[1 TimePoint]));
%                         figure(f2);plot(temp);
%                         title([num2str(bat) '-' num2str(brt)]);
%                         mask_AIF(i,j) = 0;
%                     end
                        
                end
                
                %YIJ 20170626 get rid of noisy peaks and multiple peak AIFs
                try
                %[pks2,locs2] = findpeaks(temp,'minpeakdistance',round(timepoints/20),'minpeakheight',0.6*Ctmax);
                %[pks2,locs2] = findpeaks([temp(1:bat) temp(brt:end)],'minpeakheight',0.5*Ctmax);
                %[pks2,locs2] = findpeaks(temp,'minpeakheight',0.5*Ctmax);
                %[pks2,locs2] = findpeaks([temp(1:bat) temp(brt:end)],'minpeakheight',mean(temp(temp > 0))+1.96*std(temp(temp > 0)));
                [pks2,locs2] = findpeaks([temp(1:bat) temp(brt:end)],'minpeakheight',mean(temp(ibat:brt))+2.33*std(temp(ibat:brt)));
                catch
                %disp('hi');    
                end
                
                if length(pks2) >= 1
%                     if mask_AIF(i,j)
%                         figure(1);plot(temp);
%                         figure(2);plot(reshape(smooth(images(i,j,:),5,'sgolay',3),[1 TimePoint]));
%                     end
                    mask_AIF(i,j) = 0;
                end
                     
                [pks3,locs3] = findpeaks(temp,'minpeakheight',mean(temp(ibat:brt))+2.33*std(temp(ibat:brt)));
                if isempty(pks3)
                    mask_AIF(i,j) = 0;
                end
%                 try
%                 [pks2,locs] = findpeaks(temp(bat:nmin+1),'minpeakdistance',round(length(bat:nmin+1)/10));
%                 catch
%                     disp('hi');
%                 end

                if sum(temp(bat+1:brt-1) <= 0)
                    mask_AIF(i,j) = 0;
                end
                
                cbat = temp(bat);
                if sum(temp(bat+1:brt-1) <= cbat)
                    mask_AIF(i,j) = 0;
                end

try
    indary = bat:brt;
    [maxval,maxind] = max(temp(indary));
    maxind = indary(maxind);
    [~,sortind] = sort(abs(temp(indary) - 0.5*maxval));
    sortind = indary(sortind);
    fwhm = [find(sortind < maxind,1) find(sortind > maxind,1)];
    fwhm = sortind(fwhm);
    
    tmptemp = abs(freqfilter(temp,'low',(55/200)*length(temp)));
    %[pks,locs] = findpeaks(tmptemp(fwhm(1):fwhm(2)));
    [pks,locs] = findpeaks(tmptemp,'minpeakheight',0.4*max(tmptemp(indary)));
    
%     indary = bat:brt;
%     [maxval,maxind] = max(temp(bat:brt));
%     maxind = indary(maxind);
%     window = [maxind-round((maxind-bat)*.75) maxind+round((brt-maxind)*.75)];
%     %window = [bat brt];
%     tmptemp = freqfilter(temp,'low',40);
%     [pks,locs] = findpeaks(tmptemp(window(1):window(2)));

%     figure;plot(temp(maxind-window(1):maxind+window(2)));hold on;
%     plot(locs,pks,'ro');
catch
end
                if length(pks) > 1
%                     if mask_AIF(i,j)
%                     figure(1);plot(temp);
%                        figure(2);plot(reshape(smooth(images(i,j,:),5,'sgolay',3),[1 TimePoint]));
%                        figure(3);plot(tmptemp);
%                     end
                        
                    mask_AIF(i,j)=0;
                    cond_map(i,j,s,9) = 1;
                    
                end
                
                if sum(tmptemp(bat+1:brt-1) <= 0)
                    mask_AIF(i,j) = 0;
                end
                
                if ~(length(pks) > 1) && length(pks2) >= 1
                    %figure(1);plot(temp);
                end

                cbat = tmptemp(bat);
                if sum(tmptemp(bat+1:brt-1) <= cbat)
                    mask_AIF(i,j) = 0;
                end
                
%                 
%                 try
%                 [pks2,locs2] = findpeaks(smooth(temp(bat:brt),5,'sgolay',1));
%                 catch
%                     disp('hi');
%                 end
%                 if length(pks2) > 1
%                     mask_AIF(i,j) = 0;
%                     cond_map(i,j,s,10) = 1;
%                 end
                
%                 if mask_AIF(i,j)
%                     figure(1);plot(temp);
%                        figure(2);plot(reshape(smooth(images(i,j,:),5,'sgolay',3),[1 TimePoint]));
%                        figure(3);plot(tmptemp);
%                        title([num2str(bat) '-' num2str(brt)]);
%                 end
            else
               mask_AIF(i,j)=0;
            end
        end
    end
    end
    Gd_con{s} = tmpgd_con;
    allbat{s} = tmpbat;
    filt_sig{s} = tmpfilt_sig;
    Mask_AIF(:,:,s) = mask_AIF;
end

%YIJ 20170823 Test
% Mask_AIF(138,74,1) = 1;
% Mask_AIF(139,75,1) = 1;
% Mask_AIF(151,88,1) = 1;
% Mask_AIF(157,134,1) = 1;
% Mask_AIF(158,127,1) = 1;

% % % get rid of too noisy voxels
% Mean_Std    = mean(StdsCt2Maxmap(find(Mask_AIF)));
% Stds_Std    = std(StdsCt2Maxmap(find(Mask_AIF)));
% Mask_AIF(find(StdsCt2Maxmap>Mean_Std))=0;

% max_hist    = max(StdsCt2Maxmap(find(Mask_AIF)));
% min_hist    = min(StdsCt2Maxmap(find(Mask_AIF)));
% table_Std   = min_hist:(max_hist-min_hist)/200:max_hist;
% hist_Std    = hist(StdsCt2Maxmap(find(Mask_AIF)),table_Std);

% get rid of wider bolus peak
Mean_AreaoverHight   = mean(AreaoverHightmap(find(Mask_AIF)));
Stds_AreaoverHight   = std(AreaoverHightmap(find(Mask_AIF)));
max_AreaoverHight    = max(AreaoverHightmap(find(Mask_AIF)));
min_AreaoverHight    = min(AreaoverHightmap(find(Mask_AIF)));
table_AreaoverHight  = min_AreaoverHight:(max_AreaoverHight-min_AreaoverHight)/200:max_AreaoverHight;
hist_AreaoverHight   = hist(AreaoverHightmap(find(Mask_AIF)),table_AreaoverHight);

% get rid of low Contc voxels
Mean_Conct = mean(MaxCtmap(find(Mask_AIF)));
Stds_Conct = std(MaxCtmap(find(Mask_AIF)));
max_Conct    = max(MaxCtmap(find(Mask_AIF)));
min_Conct    = 0;
table_Conct  = min_Conct:(max_Conct-min_Conct)/200:max_Conct;
hist_Conct   = hist(MaxCtmap(find(Mask_AIF)),table_Conct);

% looping for finding optimal AIF
LoopingIndex = 1;   count = 1;
HighConctIndex = 200;       LowAreaoverHightIndex = 1; 

%YIJ 20170612: This method might be better
LoopingIndex2 = 0;
n_voxel_thr = 10;
while LoopingIndex && sum(Mask_AIF(:)) > n_voxel_thr
    tempmask1 = Mask_AIF;
    tempmask2 = Mask_AIF;
    HighConctcutoff         = table_Conct(HighConctIndex);
    LowAreaoverHightcutoff  = table_AreaoverHight(LowAreaoverHightIndex);
    
    tempmask1(find(MaxCtmap<HighConctcutoff))=0;
    tempmask2(find(AreaoverHightmap>LowAreaoverHightcutoff))=0;
    
    if sum(sum(sum(tempmask1))) <= n_voxel_thr || LoopingIndex2 == 1
        HighConctIndex          = HighConctIndex -1;
    end
    if sum(sum(sum(tempmask2))) <= n_voxel_thr || LoopingIndex2 == 1
        LowAreaoverHightIndex  = LowAreaoverHightIndex + 1;
    end
    
    if sum(sum(sum(tempmask1.*tempmask2))) > n_voxel_thr
        Mask_AIF = tempmask1.*tempmask2;
        LoopingIndex = 0;
    elseif sum(sum(sum(tempmask1))) > n_voxel_thr && sum(sum(sum(tempmask2))) > n_voxel_thr
        LoopingIndex2 = 1;
    end
end

if ~sum(Mask_AIF(:))
    error('no aif found\n');
end
%YIJ 20170612: Original method commented out
% while LoopingIndex
%     tempmask = Mask_AIF;
%     HighConctcutoff         = table_Conct(HighConctIndex);
%     LowAreaoverHightcutoff  = table_AreaoverHight(LowAreaoverHightIndex);
%     
%     tempmask(find(MaxCtmap<HighConctcutoff))=0;
%     tempmask(find(AreaoverHightmap>LowAreaoverHightcutoff))=0;
%     
%     if sum(sum(sum(tempmask))) > 5
%         Mask_AIF = tempmask;
%         LoopingIndex = 0;
%     end
%     
%     HighConctIndex          = HighConctIndex -1;
%     LowAreaoverHightIndex  = LowAreaoverHightIndex + 1;
% end

%YIJ 20170823 Test
% Mask_AIF(:,:,:) = 0;
% Mask_AIF(138,74,1) = 1;
% Mask_AIF(139,75,1) = 1;
% Mask_AIF(151,88,1) = 1;
% Mask_AIF(157,134,1) = 1;
% Mask_AIF(158,127,1) = 1;

% YIJ 20170906 Implement Jessy's sampling delay fix
global sampdelayfix;
if sampdelayfix
    fprintf('Sampling delay fix\n');
    
    % YIJ: Manual AIF picking
% Mask_AIF = zeros(size(Mask_AIF));
% for s = 1:N_Slices
%     images = Images{s};
%     figure(); imshow(images(:,:,BAT+10), []); title('Make mask of brain'); mask_AIF = roipoly;
%     Mask_AIF(:,:,s) = mask_AIF;
% end
%--------------------------------------------------------------------------   
    count = 1;
for s = 1:N_Slices
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_AIF(i,j,s)
                position{count} = [i,j,s];
                count = count + 1;
            end
        end
    end
end
%--------------------------------------------------------------------------

% YIJ: Use same position from previous AIF (stored in input var, tmpdata)
if ~isempty(tmpdata)
    position = tmpdata.ROIs.positions.AIF;
    Mask_AIF = zeros(size(images,1),size(images,2),N_Slices);
    for ii = 1:length(position)
        Mask_AIF(position{ii}(1),position{ii}(2),position{ii}(3)) = 1;
    end
end
%--------------------------------------------------------------------------

    %% First: Look at slice slice with max number of AIF voxels
%Determine AIF slice position
slicePositions = zeros(1,size(position,2));
for i = 1:size(position,2)
   slicePositions(i) = position{i}(3); 
end
sliceBins = 1:N_Slices;
[N,X] = hist(slicePositions,sliceBins);
maxNindeces = find(N==max(N)); 
n_AIFslice = maxNindeces(1);

% %Adjust the AIF mask
% AIFsliceMask = zeros(size(images,1),size(images,2),N_Slices);
% AIFsliceMask(:,:,n_AIFslice) = 1;
% Mask_AIF_slice = Mask_AIF.*AIFsliceMask;
% 
% %Calculate average AIF signal and concentration
% Conct_AIF = zeros(1,TimePoint);
% Signal_AIF =  zeros(1,TimePoint);
% Signal_AIF2 = zeros(1,TimePoint);
% for i = 1:size(images,1)
%     for j = 1:size(images,2)
%         if Mask_AIF_slice(i,j,n_AIFslice)
%             %Signal_AIF = Signal_AIF + freqfilter(reshape(smooth(Images{n_AIFslice}(i,j,:),5,'sgolay',3),[1 TimePoint]),'low',10);
%             Signal_AIF = Signal_AIF + reshape(smooth(Images{n_AIFslice}(i,j,:),5,'sgolay',3),[1 TimePoint]);
%             Signal_AIF2 = Signal_AIF2 + reshape(smooth(Images{n_AIFslice}(i,j,:),.05,'sgolay',3),[1 TimePoint]);
%             Conct_AIF = Conct_AIF + reshape(Gd_con{n_AIFslice}(i,j,:),[1 TimePoint]);
%         end
%     end
% end
% Conct_AIF = Conct_AIF/sum(sum(sum(Mask_AIF_slice)));
% Conct_AIF(Conct_AIF < 0) = 0;
% Signal_AIF = Signal_AIF/sum(sum(sum(Mask_AIF_slice)));
% Signal_AIF2 = Signal_AIF2/sum(sum(sum(Mask_AIF_slice)));
% 
% % recompute AIF position
% count = 1;
% for i = 1:size(images,1)
%     for j = 1:size(images,2)
%         if Mask_AIF_slice(i,j,n_AIFslice)
%             positionAIF{count} = [i,j,n_AIFslice];
%             count = count + 1;
%         end
%     end
% end
% 
% %Determination of AIF cutoffs
% [ibat bat brt] = findBAT(Signal_AIF2);
% cutoffs_AIF = [ibat bat brt];
% 
% %% Perform gamma variate fit
% 
% %Determine true BATP and RTP
% %[bat brt] = TrueCutoffs(Conct_AIF,cutoffs_AIF);
% % shiftBATP = 0;
% % shiftRTP = -2;
% Dt = TR*1e-3;
% AIF = Conct_AIF(bat:brt);
% % AIF = Conct_AIF((bat+shiftBATP):(brt+shiftRTP));
% time = 0:Dt:(size(AIF,2)-1)*Dt;
% 
% % YIJ: take out possible clipping test
% % tmpdiff = diff(AIF);
% % upind = find(tmpdiff < 0,1);
% % downind = find(tmpdiff > 0,1,'last') + 1;
% % AIF2 = AIF([1:upind downind:end]);
% % time2 = time([1:upind downind:end]);
% 
% time2 = time;
% AIF2 = AIF;
% 
% %1) Using Grid Search algorithm to find best initial guesses for fit
% %parameters
% BestFit = gamma_variate_IDC_gridsearch(time2,AIF2);
% %Create fitted gamma curve
% SumSquares = BestFit(1);
% A = BestFit(2);
% alpha = BestFit(3);
% beta = BestFit(4);
% Fitted_AIF = gamma_variate_model_IDC(time,A,alpha,beta);
% %Plot fitted and raw AIF data points
% figure();
% plot(time,AIF,'r*',time2,AIF2,'k*',time,Fitted_AIF,'k-');
% xlabel('Time - sec');
% ylabel('Gd Concentration');
% title('AIF - Grid Search');
% legend('Average measured AIF','Gamma variate fit');
% 
% %2) Using Levenberg Marquardt least squares fitting algorithm using initial
% %guesses from the grid search step (1)
% [A,alpha,beta,RESNORM,RESIDUAL,EXITFLAG] = gamma_variate_IDC_LeastSquaresFit(time2,AIF2,A,alpha,beta);
% FitParam.A = A;
% FitParam.alpha = alpha;
% FitParam.beta = beta;
% %Create fitted gamma curve
% Fitted_AIF = gamma_variate_model_IDC(time,A,alpha,beta);
% %Plot fitted and raw AIF data points
% figure();
% plot(time,AIF,'r*',time2,AIF2,'k*',time,Fitted_AIF,'k-');
% xlabel('Time - sec');
% ylabel('Gd Concentration');
% title('AIF - Least Squares');
% legend('Average measured AIF','Gamma variate fit');
% 
% % % YIJ 20190607: remove any leading zeros
% % [~,nmax] = max(Fitted_AIF);
% % cutoffind = find(Fitted_AIF(1:nmax) > 0.0001,1);
% 
% 
% %3) Smooth AIF fit
% time_smooth = 0:Dt/N_Slices:(brt-bat)*Dt; 
% SmoothAIFfit = gamma_variate_model_IDC(time_smooth,A,alpha,beta);
% %Plot smooth, fitted and raw AIF data points
% figure();
% plot(time,AIF,'r*',time2,AIF2,'k*',time_smooth,SmoothAIFfit,'k-');
% xlabel('Time - sec');
% ylabel('Gd Concentration');
% title(['AIF - Smooth Fit: AIF slice = ' int2str(n_AIFslice)]);
% legend('Measured AIF','Gamma variate fit');
% 
% % % YIJ 20190607: remove any leading zeros
% % SmoothAIFfit = SmoothAIFfit((max([(cutoffind-1-1) 0])*N_Slices)+1:end);
% % bat = bat + cutoffind - 1 - 1;
% % time_smooth = 0:Dt/N_Slices:(brt-bat)*Dt; 
% % cutoffs_AIF = [ibat bat brt];
% 
% if RESNORM < 1e-5
%     fprintf('aif fit slice first try\n');
%     %Additional check to make sure signal has AIF shape
%     if  ( ( (max(SmoothAIFfit)-SmoothAIFfit(1)) > 0.015 ) && ( (SmoothAIFfit(length(SmoothAIFfit))-SmoothAIFfit(1)) < 0.015 ) )
%         %Return values for n_AIFslice
%         area = sum(SmoothAIFfit);
%         AIF = struct('Ct',SmoothAIFfit,'area',area,'FitParam',FitParam,'RESNORM',RESNORM,'ITP',ibat,'BATP',bat,'RTP',brt);
%         % AIF = struct('ITP',ibat,'BATP',bat,'Ct',Cmt,'RTP',brt,'area',area);
%         return
%     end
% end

%% Second: Consider all other slices if the above returns RESNORM > 1e-5,
% and pick the one with lowest RESNORM

%Determine all AIF slices
AIFSliceNumbers = X(find(N));

%%Store values from first step
% AIFsliceMasks{find(AIFSliceNumbers==n_AIFslice)} = Mask_AIF_slice;
% AIFpositions{find(AIFSliceNumbers==n_AIFslice)} = positionAIF;
% AIFsignals{find(AIFSliceNumbers==n_AIFslice)} = Signal_AIF;
% AIFconcts{find(AIFSliceNumbers==n_AIFslice)} = Conct_AIF;
% AIFcutoffs{find(AIFSliceNumbers==n_AIFslice)} = cutoffs_AIF;
% AIFsmoothFits{find(AIFSliceNumbers==n_AIFslice)} = SmoothAIFfit;
% AIFfitParams{find(AIFSliceNumbers==n_AIFslice)} = FitParam;
% % YIJ 20180910 change initial resnorm to Inf
% %Resnorms = zeros(1,size(AIFSliceNumbers,2));
% Resnorms = ones(1,size(AIFSliceNumbers,2)).*Inf;
% Resnorms(find(AIFSliceNumbers==n_AIFslice)) = RESNORM;
% ibats = zeros(1,size(AIFSliceNumbers,2));
% ibats(find(AIFSliceNumbers==n_AIFslice)) = ibat;
% bats = zeros(1,size(AIFSliceNumbers,2));
% bats(find(AIFSliceNumbers==n_AIFslice)) = bat;
% brts = zeros(1,size(AIFSliceNumbers,2));
% brts(find(AIFSliceNumbers==n_AIFslice)) = brt;
n_AIFslice_max = n_AIFslice;
% tmptimes{find(AIFSliceNumbers==n_AIFslice)} = time;
% tmptimes2{find(AIFSliceNumbers==n_AIFslice)} = time2;
% tmpaifs{find(AIFSliceNumbers==n_AIFslice)} = AIF;
% tmpaifs2{find(AIFSliceNumbers==n_AIFslice)} = AIF2;
% tmpsmoothtimes{find(AIFSliceNumbers==n_AIFslice)} = time_smooth;
% tmpsmoothaifs{find(AIFSliceNumbers==n_AIFslice)} = SmoothAIFfit;

%Looping over all remaining AIF slices
for s = 1:size(AIFSliceNumbers,2)
    n_AIFslice = AIFSliceNumbers(s);
    %Skip n_AIFslice_max from step 1
%     if n_AIFslice == n_AIFslice_max
%         continue
%     end
    %Adjust the AIF mask
    AIFsliceMask = zeros(size(images,1),size(images,2),N_Slices);
    AIFsliceMask(:,:,n_AIFslice) = 1;
    Mask_AIF_slice = Mask_AIF.*AIFsliceMask;
    %Calculate average AIF signal and concentration
    Conct_AIF = zeros(1,TimePoint);
    Signal_AIF =  zeros(1,TimePoint);
    Signal_AIF2 =  zeros(1,TimePoint);
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_AIF_slice(i,j,n_AIFslice)
                %Signal_AIF = Signal_AIF + freqfilter(reshape(smooth(Images{n_AIFslice}(i,j,:),5,'sgolay',3),[1 TimePoint]),'low',10);;
                %Signal_AIF = Signal_AIF + reshape(smooth(Images{n_AIFslice}(i,j,:),5,'sgolay',3),[1 TimePoint]);
                %Signal_AIF2 = Signal_AIF2 + reshape(smooth(Images{n_AIFslice}(i,j,:),max(sqrt(1/(TR/1000))*5,5),'sgolay',3),[1 TimePoint]);
                Conct_AIF = Conct_AIF + reshape(Gd_con{n_AIFslice}(i,j,:),[1 TimePoint]);
                %figure;plot(Signal_AIF);
                
                % 20210422 YIJ: remove smoothing because of arrival time shift
                Signal_AIF = Signal_AIF + reshape(Images{n_AIFslice}(i,j,:),[1 timepoints]);
                % 20210329 YIJ: do spline smoothing
%                 tmpfit = fit([1:timepoints]',squeeze(Images{n_AIFslice}(i,j,:)),'smoothingspline','smoothingparam',.8);
%                 Signal_AIF = Signal_AIF + reshape(tmpfit(1:timepoints),[1 timepoints]);
%                 Signal_AIF(Signal_AIF < 0) = 0;
                Signal_AIF2 = Signal_AIF;

            end
        end
    end
    Conct_AIF = Conct_AIF/sum(sum(sum(Mask_AIF_slice)));
    Conct_AIF(Conct_AIF < 0) = 0;
    Signal_AIF = Signal_AIF/sum(sum(sum(Mask_AIF_slice)));
    Signal_AIF2 = Signal_AIF2/sum(sum(sum(Mask_AIF_slice)));

    %Determination of AIF cutoffs
    [ibat bat brt] = findBAT_test(Signal_AIF2);
    
    % YIJ 20180910
    if bat < 1 || brt >= 199
        figure;plot(Signal_AIF2);title(['bad aif on slc ' num2str(n_AIFslice)]);
        continue;
    end
    cutoffs_AIF = [ibat bat brt];
    
    % recompute AIF position 
    count = 1;
    clear positionAIF
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_AIF_slice(i,j,n_AIFslice)
                positionAIF{count} = [i,j,n_AIFslice];
                count = count + 1;
            end
        end
    end
    
    %%%Perform gamma variate fit
    
    %Determine true BATP and RTP
    %[bat brt] = TrueCutoffs(Conct_AIF,cutoffs_AIF);
%     shiftBATP = 0;
%     shiftRTP = -2;
    Dt = TR*1e-3;
    AIF = Conct_AIF(bat:brt);
%     AIF = Conct_AIF((bat+shiftBATP):(brt+shiftRTP));
    time = 0:Dt:(size(AIF,2)-1)*Dt;

    % take out possible clipping test
%     tmpdiff = diff(AIF);
%     upind = find(tmpdiff < 0,1);
%     downind = find(tmpdiff > 0,1,'last') + 1;
%     AIF2 = AIF([1:upind downind:end]);
%     time2 = time([1:upind downind:end]);
    
    time2 = time;
    AIF2 = AIF;
    
    %1) Using Grid Search algorithm to find best initial guesses for fit
    %parameters
    BestFit = gamma_variate_IDC_gridsearch(time2,AIF2);
    %Create fitted gamma curve
    SumSquares = BestFit(1);
    A = BestFit(2);
    alpha = BestFit(3);
    beta = BestFit(4);
    Fitted_AIF = gamma_variate_model_IDC(time,A,alpha,beta);
    %Plot fitted and raw AIF data points
    if flg_plot
    figure();
    plot(time,AIF,'r*',time2,AIF2,'k*',time,Fitted_AIF,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('AIF - Grid Search');
    legend('Average measured AIF','Gamma variate fit');
    end

    %2) Using Levenberg Marquardt least squares fitting algorithm using initial
    %guesses from the grid search step (1)
    [A,alpha,beta,RESNORM,RESIDUAL,EXITFLAG] = gamma_variate_IDC_LeastSquaresFit(time2,AIF2,A,alpha,beta);
    FitParam.A = A;
    FitParam.alpha = alpha;
    FitParam.beta = beta;
    %Create fitted gamma curve
    Fitted_AIF = gamma_variate_model_IDC(time,A,alpha,beta);
    %Plot fitted and raw AIF data points
    if flg_plot
    figure();
    plot(time,AIF,'r*',time2,AIF2,'k*',time,Fitted_AIF,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('AIF - Least Squares');
    legend('Average measured AIF','Gamma variate fit');
    end
    
% % YIJ 20190607: remove any leading zeros
% [~,nmax] = max(Fitted_AIF);
% cutoffind = find(Fitted_AIF(1:nmax) > 0.0001,1);
    
    % Smooth AIF fit
    time_smooth = 0:Dt/N_Slices:(brt-bat)*Dt;
    SmoothAIFfit = gamma_variate_model_IDC(time_smooth,A,alpha,beta);
    %Plot smooth, fitted and raw AIF data points
    if flg_plot
    figure();
    plot(time,AIF,'r*',time2,AIF2,'k*',time_smooth,SmoothAIFfit,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('AIF - Smooth Fit');
    legend('Measured AIF','Gamma variate fit');
    end
    
%     % YIJ 20190607: remove any leading zeros
% SmoothAIFfit = SmoothAIFfit((max([(cutoffind-1-1) 0])*N_slices)+1:end);
% bat = bat + cutoffind - 1 - 1;
% time_smooth = 0:Dt/N_Slices:(brt-bat)*Dt;
% cutoffs_AIF = [ibat bat brt];

    AIFsliceMasks{find(AIFSliceNumbers==n_AIFslice)} = Mask_AIF_slice;
    AIFpositions{find(AIFSliceNumbers==n_AIFslice)} = positionAIF;
    AIFsignals{find(AIFSliceNumbers==n_AIFslice)} = Signal_AIF;
    AIFconcts{find(AIFSliceNumbers==n_AIFslice)} = Conct_AIF;
    AIFcutoffs{find(AIFSliceNumbers==n_AIFslice)} = cutoffs_AIF;
    AIFsmoothFits{find(AIFSliceNumbers==n_AIFslice)} = SmoothAIFfit;
    AIFfitParams{find(AIFSliceNumbers==n_AIFslice)} = FitParam;
    Resnorms(find(AIFSliceNumbers==n_AIFslice)) = RESNORM;
    ibats(find(AIFSliceNumbers==n_AIFslice)) = ibat;
    bats(find(AIFSliceNumbers==n_AIFslice)) = bat;
    brts(find(AIFSliceNumbers==n_AIFslice)) = brt;
    tmptimes{find(AIFSliceNumbers==n_AIFslice)} = time;
tmptimes2{find(AIFSliceNumbers==n_AIFslice)} = time2;
tmpaifs{find(AIFSliceNumbers==n_AIFslice)} = AIF;
tmpaifs2{find(AIFSliceNumbers==n_AIFslice)} = AIF2;
tmpsmoothtimes{find(AIFSliceNumbers==n_AIFslice)} = time_smooth;
tmpsmoothaifs{find(AIFSliceNumbers==n_AIFslice)} = SmoothAIFfit;
        
%     if RESNORM < 1e-5
%         break
%     end

end

% %If RESNORM does not go below 1e-5, search for slice with lowest RESNORM
% if RESNORM >= 1e-5
%     n_AIFslice = AIFSliceNumbers(find(Resnorms==min(Resnorms)));
% end

%Loop to determine n_AIFslice
Resnorms_temp = Resnorms;
LoopingIndex = 1;
while LoopingIndex
    n_AIFslice = AIFSliceNumbers(find(Resnorms_temp==min(Resnorms_temp)));
    SmoothAIFfit = AIFsmoothFits{find(AIFSliceNumbers==n_AIFslice)};
    %Check to make sure n_AIFslice gives an actual AIF signal in shape
    if ( ( (max(SmoothAIFfit)-SmoothAIFfit(1)) > 0.015 ) && ( (SmoothAIFfit(length(SmoothAIFfit))-SmoothAIFfit(1)) < 0.015 ) )
        break;
    end  
    Resnorms_temp(find(AIFSliceNumbers==n_AIFslice)) = 999999;
    %Exit point
    if (min(Resnorms_temp) == 999999)
        n_AIFslice = n_AIFslice_max;
        break;
    end
end

if ~isempty(tmpdata)
    n_AIFslice = tmpdata.ROIs.positions.n_slice_AIF;
end
%n_AIFslice = maxNindeces(1);

%n_AIFslice = 4;
%%Return values for n_AIFslice
Mask_AIF_slice = AIFsliceMasks{find(AIFSliceNumbers==n_AIFslice)};
positionAIF = AIFpositions{find(AIFSliceNumbers==n_AIFslice)};
Signal_AIF = AIFsignals{find(AIFSliceNumbers==n_AIFslice)};
Conct_AIF = AIFconcts{find(AIFSliceNumbers==n_AIFslice)};
cutoffs_AIF = AIFcutoffs{find(AIFSliceNumbers==n_AIFslice)};

FitParam = AIFfitParams{find(AIFSliceNumbers==n_AIFslice)};
RESNORM = Resnorms(find(AIFSliceNumbers==n_AIFslice));
ibat = ibats(find(AIFSliceNumbers==n_AIFslice));
bat = bats(find(AIFSliceNumbers==n_AIFslice));
brt = brts(find(AIFSliceNumbers==n_AIFslice));
SmoothAIFfit = AIFsmoothFits{find(AIFSliceNumbers==n_AIFslice)};

tmptime = tmptimes{find(AIFSliceNumbers==n_AIFslice)};
tmptime2 = tmptimes2{find(AIFSliceNumbers==n_AIFslice)};
tmpaif = tmpaifs{find(AIFSliceNumbers==n_AIFslice)};
tmpaif2 = tmpaifs2{find(AIFSliceNumbers==n_AIFslice)};
tmpsmoothtime = tmpsmoothtimes{find(AIFSliceNumbers==n_AIFslice)};
tmpsmoothaif = tmpsmoothaifs{find(AIFSliceNumbers==n_AIFslice)};

figure;plot(Signal_AIF);
    title(['aif: ' num2str(bat) '-' num2str(brt) ' TR: ' num2str(TR/1000)]);
    %exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_aif(auto).jpg'], 'resolution', 300);
    
figure();
plot(tmptime,tmpaif,'r*',tmptime2,tmpaif2,'k*',tmpsmoothtime,tmpsmoothaif,'k-');
xlabel('Time - sec');
ylabel('Gd Concentration');
title(['AIF - Smooth Fit - chosen slice ' num2str(n_AIFslice)]);
legend('Measured AIF','Gamma variate fit');
%exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_aifFit(auto).jpg'], 'resolution', 300);
pause(0.1);

%yv(Images{n_AIFslice}(:,:,1), [], 'ol', Mask_AIF_slice(:,:,n_AIFslice))
figure;imshow(Images{n_AIFslice}(:,:,1),[],'initialmagnification','fit');
    roir = ones(size(Images{n_AIFslice},[1 2]));
    roib = zeros(size(Images{n_AIFslice},[1 2]));
    roig = zeros(size(Images{n_AIFslice},[1 2]));
    roic = cat(3, roir, roig, roib);
    hold on; imh = imshow(roic);
    set(imh, 'alphadata', 0.5*Mask_AIF_slice(:,:,n_AIFslice));

%exportgraphics(gca, ['..\images\autovein debug\' currtime '_' caseid '_aifloc(auto).jpg'], 'resolution', 300);

area = sum(SmoothAIFfit);
AIF = struct('Ct',SmoothAIFfit,'area',area,'FitParam',FitParam,'RESNORM',RESNORM,'ITP',ibat,'BATP',bat,'RTP',brt,'signal',Signal_AIF,'conc',Conct_AIF);
% AIF = struct('ITP',ibat,'BATP',bat,'Ct',Cmt,'RTP',brt,'area',area);

else
% accumulate
Conct_AIF = zeros(1,TimePoint);
Signal_AIF =  zeros(1,TimePoint);
Signal_AIF2 =  zeros(1,TimePoint);
for s = 1:N_Slices
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_AIF(i,j,s)
                %YIJ 20170614 smooth\
                %Signal_AIF = Signal_AIF + freqfilter(reshape(smooth(Images{s}(i,j,:),5,'sgolay',3),[1 TimePoint]),'low',10);
                %Signal_AIF = Signal_AIF + reshape(smooth(Images{s}(i,j,:),5,'sgolay',3),[1 TimePoint]);
                %Signal_AIF2 = Signal_AIF2 + reshape(smooth(Images{s}(i,j,:),max(sqrt(1/(TR/1000))*5,5),'sgolay',3),[1 TimePoint]);
                %Signal_AIF = Signal_AIF + reshape(Images{s}(i,j,:),[1 TimePoint]);
                Conct_AIF = Conct_AIF + reshape(Gd_con{s}(i,j,:),[1 TimePoint]);
                
                % 20210422 YIJ: remove smoothing because of arrival time shift
                Signal_AIF = Signal_AIF + reshape(Images{s}(i,j,:),[1 timepoints]);
                % 20210329 YIJ: do spline smoothing
%                 tmpfit = fit([1:timepoints]',Images{s}(i,j,:)','smoothingspline','smoothingparam',.8);
%                 Signal_AIF = Signal_AIF + reshape(tmpfit(1:timepoints),[1 timepoints]);
%                 Signal_AIF(Signal_AIF < 0) = 0;
                Signal_AIF2 = Signal_AIF;
            end
        end
    end
end
Conct_AIF = Conct_AIF/sum(sum(sum(Mask_AIF)));
Conct_AIF(Conct_AIF < 0) = 0;
Signal_AIF = Signal_AIF/sum(sum(sum(Mask_AIF)));
Signal_AIF2 = Signal_AIF2/sum(sum(sum(Mask_AIF)));

% position
count = 1;
for s = 1:N_Slices
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_AIF(i,j,s)
                position{count} = [i,j,s];
                count = count + 1;
            end
        end
    end
end

%Determination of AIF cutoffs
[ibat bat brt] = findBAT_test(Signal_AIF2);
cutoffs_AIF = [ibat bat brt];

%YIJ 20170620
figure;plot(Signal_AIF);title([num2str(ibat) ' ' num2str(bat) ' ' num2str(brt)]);
end

end
    
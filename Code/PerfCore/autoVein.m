function [Vein,Mask_Vein_slice,positionVein,n_Veinslice] = autoVein(path_DSC,N_meas,N_Slices,AIFslice)
%autoVein Automatically finds a vein
%   Finds a vein similar to autoAIF
%
% Author: Yong Ik Jeong
% Date: 2020-05-15
% Changelog:
%   - 20200515 YIJ: Initial version copied from autoaif_wy_Philips_v2
%   - 20200612 YIJ: changed 2_1 and 2_2 to +2 from +1 for late arrival
%   - 20200616 YIJ: changed 2_1 and 2_2 to +3

flg_plot = 0;

global caseid;
global currtime;
% if ~exist('caseid','var')
%     caseid = split(path_DSC,{'\','/'});
%     caseid = caseid{end-3};
% end
% currtime = datestr(now,'yyyymmddHHMMSS');

header = dicominfo([path_DSC '\1.dcm'],'Dictionary','dicom-dict.txt'); %Upon the upgrade of the dicom dictionary

for  s = 1:N_Slices
    tempimage = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(s-1)*N_meas),(s*N_meas),1))); 
        Images{s} = tempimage;
end

% 20210326 YIJ: spatially smooth DSC (xy and z)
Images = smoothDSC(Images);

TR = header.RepetitionTime; %ms    
TE = header.EchoTime;       %ms
timepoints = size(Images{1},3);


NmaskNS = 0;
batlist = [];
brtlist = [];
nminlist = [];

for s = 1:N_Slices

    images = Images{s};

    mask_s = automaskns(imfilter(images(:,:,1),fspecial('gaussian',[3 3],2)));
      
    sum_signal = zeros(1,size(images,3));
    for i = 1:size(images,1)
    for j = 1:size(images,2)
    if mask_s(i,j)
       
        %tmpary = reshape(smooth(images(i,j,:),max(sqrt(1/(TR/1000))*5,5),'sgolay',3),[1 size(images,3)]);
        
        % 20210422 YIJ: remove smoothing because of arrival time shift
        tmpary = reshape(images(i,j,:),[1 timepoints]);
        % 20210329 YIJ: do spline smoothing
%         tmpfit = fit([1:timepoints]',images(i,j,:)','smoothingspline','smoothingparam',.8);
%         tmpary = reshape(tmpfit(1:timepoints),[1 timepoints]);
%         tmpary(tmpary < 0) = 0;
        
        sum_signal = sum_signal+ tmpary;
        
        % 20210210 YIJ: collect list of bat, brt, and nmin
        [~,tmpbat,tmpbrt] = findBAT_test(tmpary);
        [~,tmpnmin] = min(tmpary);
        batlist = [batlist tmpbat];
        brtlist = [brtlist tmpbrt];
        nminlist = [nminlist tmpnmin];
    end
    end
    end

    Sum_Signal(s,:) =  sum_signal;
    NmaskNS = NmaskNS + sum(sum(mask_s));
end

Mean_Signal = sum(Sum_Signal)/NmaskNS;

time = 0:TR*1e-3:(size(images,3)-1)*TR*1e-3;

% pre definition
TimePoint           = size(images,3);
MaxCtmap            = zeros(size(images,1),size(images,2),N_Slices);
% IntegralCt2BRTmap   = zeros(size(images,1),size(images,2),N_Slices);
AreaoverHightmap    = zeros(size(images,1),size(images,2),N_Slices);
% StdsCt2Maxmap       = zeros(size(images,1),size(images,2),N_Slices);
Mask_Vein            = zeros(size(images,1),size(images,2),N_Slices);

cond_map = zeros(size(images,1),size(images,2),N_Slices,10);
tmpneg_map = zeros(size(images,1),size(images,2),N_Slices);

% 20210223 YIJ: Use AIF instead of global average as per Tim's suggestion
Mean_Signal = AIFslice.signal;

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

for s = 1:N_Slices
    % read images
    images = Images{s};

    %Create new mask 
    %YIJ 20170614 smooth mask
    mask_Vein = automaskns(imfilter(images(:,:,1),fspecial('gaussian',[3 3],2)));
    
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
        %if i == 117 && j == 63 && s == 26
        if i == 37 && j == 41 && s == 5
            %disp('debug');
        end
        if mask_Vein(i,j)

            %YIJ 20170614 smooth
            %tempsig = reshape(smooth(images(i,j,:),5,'sgolay',3),[1 TimePoint]);
            %tempsig2 = reshape(smooth(images(i,j,:),max(sqrt(1/(TR/1000))*5,5),'sgolay',3),[1 TimePoint]);            
            
            % 20210422 YIJ: remove smoothing because of arrival time shift
            tempsig = reshape(images(i,j,:),[1 timepoints]);
            % 20210329 YIJ: do spline smoothing
%             tmpfit = fit([1:timepoints]',images(i,j,:)','smoothingspline','smoothingparam',.8);
%             tempsig = reshape(tmpfit(1:timepoints),[1 timepoints]);
            tempsig2 = tempsig;
                        
            %temp2 = temp;
            %YIJ 20170626 smoothing can cause negative numbers sometimes
            tempsig(tempsig < 0) = 0;
            tempsig2(tempsig2 < 0) = 0;           
            
            meanS0 = mean(tempsig(iBAT:BAT));
            stdS0 = std(tempsig(iBAT:BAT));
            
%             if meanS0 < 400
%                 mask_AIF(i,j) = 0;
%             end
               
            [smax nmax] = max(tempsig);   % temp = signal
            [smin nmin] = min(tempsig);   % temp = signal

            % calcualte iBAT, BAT & BRT
            [ibat bat brt] = findBAT_test(tempsig2);
            
            if ibat & bat & brt
            
                % Low baseline signal is probably wrong
                if mean(tempsig(ibat:bat)) < 200
                    mask_Vein(i,j) = 0;
                end
                
                % make Gd_concent maps
                temp = -1/TE*log(tempsig./mean(tempsig(ibat:bat))); % temp = conct.
                temp(find(~isfinite(temp)))=0;
                tmpgd_con(i,j,:) = temp;
                %YIJ 20170608
                temp(temp < 0) = 0;
                
                temp2 = -1/TE*log(tempsig2./mean(tempsig2(ibat:bat))); % temp = conct.
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
                    %20210204 YIJ comment out 
                    %mask_Vein(i,j)=0;
                    cond_map(i,j,s,1) = 1;
                end
                % 2_1. remove early BAT based on Maximum concentration
                if nmin < Nmin + 3
                    mask_Vein(i,j)=0;
                    cond_map(i,j,s,2) = 1;
                end
                % 2_2. remove early BAT based on BAT
                if bat < BAT + 3
                    mask_Vein(i,j)=0;
                    cond_map(i,j,s,3) = 1;
                end
                if bat < AIFslice.BATP + 3
                    mask_Vein(i,j)=0;
                end
                % 3. remove early recirculation
                if brt < BRT + 2
                    mask_Vein(i,j)=0;
                    cond_map(i,j,s,4) = 1;
                end
                % 4. remove smaller integral of Conct. up to recir.
                if sum(temp(bat:brt)) < sum(Mean_Conct(BAT:BRT))
                    mask_Vein(i,j)=0;
                    cond_map(i,j,s,5) = 1;
                end
                % 5. remove lower highest peak.
                if Ctmax < max(Mean_Conct)
                    mask_Vein(i,j)=0;
                    cond_map(i,j,s,6) = 1;
                end
                % 6. remove narrow bolus.
                if sum(temp(bat:brt))/Ctmax < sum(Mean_Conct(BAT:BRT))/MeanCtmax
                    mask_Vein(i,j)=0;
                    cond_map(i,j,s,7) = 1;
                end
                
                %YIJ 20170818 bolus peak should be between bat and brt
                if nmin < bat || nmin > brt
                    mask_Vein(i,j) = 0;
                    cond_map(i,j,s,8) = 1;                        
                end
                
                % 20210204 YIJ: disregard very late bolus
                % can't accurately sample the bolus AND it's usually not a
                % vein
                if nmin > TimePoint * 0.8
                    mask_Vein(i,j) = 0;
                    cond_map(i,j,s,9) = 1;
                end
                
                % 20210210 YIJ: get rid of bat less than 90th percentile
                if bat < prctile(batlist, 80)
                    mask_Vein(i,j) = 0;
                    cond_map(i,j,s,10) = 1;
                end
                
                % 20210210 YIJ: try remove leakage signals
                if mean(tempsig(ibat:ibat+5)) < 0.8 * mean(tempsig(end-6:end-1))
                    mask_Vein(i,j) = 0;
                    cond_map(i,j,s,11) = 1;
                end
                
                %YIJ 20170626 get rid of noisy peaks and multiple peak AIFs
                try
                [pks2,locs2] = findpeaks([temp(1:bat) temp(brt:end)],'minpeakheight',mean(temp(ibat:brt))+2.33*std(temp(ibat:brt)));
                catch
                %disp('hi');    
                end
                
                if length(pks2) >= 1
                    mask_Vein(i,j) = 0;
                end
                     
                [pks3,locs3] = findpeaks(temp,'minpeakheight',mean(temp(ibat:brt))+2.33*std(temp(ibat:brt)));
                if isempty(pks3)
                    mask_Vein(i,j) = 0;
                end

                if sum(temp(bat+1:brt-1) <= 0)
                    mask_Vein(i,j) = 0;
                end
                
                cbat = temp(bat);
                if sum(temp(bat+1:brt-1) <= cbat)
                    mask_Vein(i,j) = 0;
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

catch
end
                if length(pks) > 1                        
                    mask_Vein(i,j)=0;
                    cond_map(i,j,s,12) = 1;                    
                end
                
                if sum(tmptemp(bat+1:brt-1) <= 0)
                    mask_Vein(i,j) = 0;
                end
                
                if ~(length(pks) > 1) && length(pks2) >= 1
                    %figure(1);plot(temp);
                end

                cbat = tmptemp(bat);
                if sum(tmptemp(bat+1:brt-1) <= cbat)
                    mask_Vein(i,j) = 0;
                end
                
%                 if isempty(pks) || isempty(pks2) || isempty(pks3)
%                     mask_AIF(i,j) = 0;
%                 end
            else
               mask_Vein(i,j)=0;
            end
        end
    end
    end
    Gd_con{s} = tmpgd_con;
    allbat{s} = tmpbat;
    filt_sig{s} = tmpfilt_sig;
    Mask_Vein(:,:,s) = mask_Vein;
end
% tmpind = 1:N_Slices;
% tmpind(aifslc) = [];
% Mask_Vein(:,:,tmpind) = 0;

% get rid of wider bolus peak
Mean_AreaoverHight   = mean(AreaoverHightmap(find(Mask_Vein)));
Stds_AreaoverHight   = std(AreaoverHightmap(find(Mask_Vein)));
max_AreaoverHight    = max(AreaoverHightmap(find(Mask_Vein)));
min_AreaoverHight    = min(AreaoverHightmap(find(Mask_Vein)));
table_AreaoverHight  = min_AreaoverHight:(max_AreaoverHight-min_AreaoverHight)/200:max_AreaoverHight;
hist_AreaoverHight   = hist(AreaoverHightmap(find(Mask_Vein)),table_AreaoverHight);

% get rid of low Contc voxels
Mean_Conct = mean(MaxCtmap(find(Mask_Vein)));
Stds_Conct = std(MaxCtmap(find(Mask_Vein)));
max_Conct    = max(MaxCtmap(find(Mask_Vein)));
min_Conct    = 0;
table_Conct  = min_Conct:(max_Conct-min_Conct)/200:max_Conct;
hist_Conct   = hist(MaxCtmap(find(Mask_Vein)),table_Conct);

% looping for finding optimal AIF
LoopingIndex = 1;   count = 1;
HighConctIndex = 200;       HighAreaoverHightIndex = 200; 

%YIJ 20170612: This method might be better
LoopingIndex2 = 0;
n_voxel_thr = 5;
while LoopingIndex && sum(Mask_Vein(:)) > n_voxel_thr
    tempmask1 = Mask_Vein;
    tempmask2 = Mask_Vein;
    HighConctcutoff         = table_Conct(HighConctIndex);
    HighAreaoverHightcutoff  = table_AreaoverHight(HighAreaoverHightIndex);
    
    % 20200219 Try change this
    tempmask1(find(MaxCtmap<HighConctcutoff))=0;
    tempmask2(find(AreaoverHightmap<HighAreaoverHightcutoff))=0;
    
    if sum(sum(sum(tempmask1))) <= n_voxel_thr || LoopingIndex2 == 1
        HighConctIndex          = HighConctIndex -1;
    end
    if sum(sum(sum(tempmask2))) <= n_voxel_thr || LoopingIndex2 == 1
        HighAreaoverHightIndex  = HighAreaoverHightIndex - 1;
    end
    
    if sum(sum(sum(tempmask1.*tempmask2))) > n_voxel_thr
        Mask_Vein = tempmask1.*tempmask2;
        LoopingIndex = 0;
    elseif sum(sum(sum(tempmask1))) > n_voxel_thr && sum(sum(sum(tempmask2))) > n_voxel_thr
        LoopingIndex2 = 1;
    end
end

if ~sum(Mask_Vein(:))
    error('no vein found\n');
end

% debug
dsc = cat(4, Images{:});
plotMapSignal(dsc, Mask_Vein);

% Find slice with most chosen voxels
    count = 1;
for s = 1:N_Slices
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_Vein(i,j,s)
                position{count} = [i,j,s];
                count = count + 1;
            end
        end
    end
end

slicePositions = zeros(1,size(position,2));
for i = 1:size(position,2)
   slicePositions(i) = position{i}(3); 
end
sliceBins = 1:N_Slices;
[N,X] = hist(slicePositions,sliceBins);
maxNindeces = find(N==max(N)); 
n_Veinslice = maxNindeces(1);

VeinSliceNumbers = X(find(N));
n_Veinslice_max = n_Veinslice;

% Loop through chosen slices to find lowest RESNORM
for s = 1:size(VeinSliceNumbers,2)
    n_Veinslice = VeinSliceNumbers(s);
    %Skip n_Veinslice_max from step 1
%     if n_Veinslice == n_Veinslice_max
%         continue
%     end
    %Adjust the Vein mask
    VeinsliceMask = zeros(size(images,1),size(images,2),N_Slices);
    VeinsliceMask(:,:,n_Veinslice) = 1;
    Mask_Vein_slice = Mask_Vein.*VeinsliceMask;
    %Calculate average Vein signal and concentration
    Conct_Vein = zeros(1,TimePoint);
    Signal_Vein =  zeros(1,TimePoint);
    Signal_Vein2 =  zeros(1,TimePoint);
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_Vein_slice(i,j,n_Veinslice)
                %Signal_Vein = Signal_Vein + freqfilter(reshape(smooth(Images{n_Veinslice}(i,j,:),5,'sgolay',3),[1 TimePoint]),'low',10);;
                %Signal_Vein = Signal_Vein + reshape(smooth(Images{n_Veinslice}(i,j,:),5,'sgolay',3),[1 TimePoint]);
                %Signal_Vein2 = Signal_Vein2 + reshape(smooth(Images{n_Veinslice}(i,j,:),max(sqrt(1/(TR/1000))*5,5),'sgolay',3),[1 TimePoint]);
                Conct_Vein = Conct_Vein + reshape(Gd_con{n_Veinslice}(i,j,:),[1 TimePoint]);
                %figure;plot(Signal_Vein);
                
                % 20210422 YIJ: remove smoothing because of arrival time shift
                Signal_Vein = Signal_Vein + reshape(Images{n_Veinslice}(i,j,:),[1 timepoints]);
                % 20210329 YIJ: do spline smoothing
%                 tmpfit = fit([1:timepoints]',Images{n_AIFslice}(i,j,:)','smoothingspline','smoothingparam',.8);
%                 Signal_Vein = Signal_Vein + reshape(tmpfit(1:timepoints),[1 timepoints]);
%                 Signal_Vein(Signal_Vein < 0) = 0;
                Signal_Vein2 = Signal_Vein;

            end
        end
    end
    Conct_Vein = Conct_Vein/sum(sum(sum(Mask_Vein_slice)));
    Conct_Vein(Conct_Vein < 0) = 0;
    Signal_Vein = Signal_Vein/sum(sum(sum(Mask_Vein_slice)));
    Signal_Vein2 = Signal_Vein2/sum(sum(sum(Mask_Vein_slice)));
    
[ibat, bat, brt] = findBAT_test(Signal_Vein2);
% [tmpibat,tmpbat,tmpbrt] = findBAT(Signal_Vein);
% if abs(bat-tmpbat)*Dt > 1
%     if tmpbat > bat
%         ibat = tmpibat;
%         bat = tmpbat;
%         brt = tmpbrt;
%     end
%     
% end
cutoffs_Vein = [ibat bat brt];

% recompute Vein position 
    count = 1;
    clear positionVein
    for i = 1:size(images,1)
        for j = 1:size(images,2)
            if Mask_Vein_slice(i,j,n_Veinslice)
                positionVein{count} = [i,j,n_Veinslice];
                count = count + 1;
            end
        end
    end
    
    %%%Perform gamma variate fit
    
    %Determine true BATP and RTP
    %[bat brt] = TrueCutoffs(Conct_Vein,cutoffs_Vein);
%     shiftBATP = 0;
%     shiftRTP = -2;
    Dt = TR*1e-3;
    Vein = Conct_Vein(bat:brt);
%     Vein = Conct_Vein((bat+shiftBATP):(brt+shiftRTP));
    time = 0:Dt:(size(Vein,2)-1)*Dt;

    % take out possible clipping test
%     tmpdiff = diff(Vein);
%     upind = find(tmpdiff < 0,1);
%     downind = find(tmpdiff > 0,1,'last') + 1;
%     Vein2 = Vein([1:upind downind:end]);
%     time2 = time([1:upind downind:end]);
    
    time2 = time;
    Vein2 = Vein;
    
    %1) Using Grid Search algorithm to find best initial guesses for fit
    %parameters
    BestFit = gamma_variate_IDC_gridsearch(time2,Vein2);
    %Create fitted gamma curve
    SumSquares = BestFit(1);
    A = BestFit(2);
    alpha = BestFit(3);
    beta = BestFit(4);
    Fitted_Vein = gamma_variate_model_IDC(time,A,alpha,beta);
    %Plot fitted and raw Vein data points
    if flg_plot
    figure();
    plot(time,Vein,'r*',time2,Vein2,'k*',time,Fitted_Vein,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('Vein - Grid Search');
    legend('Average measured Vein','Gamma variate fit');
    end

    %2) Using Levenberg Marquardt least squares fitting algorithm using initial
    %guesses from the grid search step (1)
    [A,alpha,beta,RESNORM,RESIDUAL,EXITFLAG] = gamma_variate_IDC_LeastSquaresFit(time2,Vein2,A,alpha,beta);
    FitParam.A = A;
    FitParam.alpha = alpha;
    FitParam.beta = beta;
    %Create fitted gamma curve
    Fitted_Vein = gamma_variate_model_IDC(time,A,alpha,beta);
    %Plot fitted and raw Vein data points
    if flg_plot
    figure();
    plot(time,Vein,'r*',time2,Vein2,'k*',time,Fitted_Vein,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('Vein - Least Squares');
    legend('Average measured Vein','Gamma variate fit');
    end
    
% % YIJ 20190607: remove any leading zeros
% [~,nmax] = max(Fitted_Vein);
% cutoffind = find(Fitted_Vein(1:nmax) > 0.0001,1);
    
    % Smooth Vein fit
    time_smooth = 0:Dt/N_Slices:(brt-bat)*Dt;
    SmoothVeinfit = gamma_variate_model_IDC(time_smooth,A,alpha,beta);
    %Plot smooth, fitted and raw Vein data points
    if flg_plot
    figure();
    plot(time,Vein,'r*',time2,Vein2,'k*',time_smooth,SmoothVeinfit,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('Vein - Smooth Fit');
    legend('Measured Vein','Gamma variate fit');
    end
    
%     % YIJ 20190607: remove any leading zeros
% SmoothVeinfit = SmoothVeinfit((max([(cutoffind-1-1) 0])*N_slices)+1:end);
% bat = bat + cutoffind - 1 - 1;
% time_smooth = 0:Dt/N_Slices:(brt-bat)*Dt;
% cutoffs_Vein = [ibat bat brt];

    VeinsliceMasks{find(VeinSliceNumbers==n_Veinslice)} = Mask_Vein_slice;
    Veinpositions{find(VeinSliceNumbers==n_Veinslice)} = positionVein;
    Veinsignals{find(VeinSliceNumbers==n_Veinslice)} = Signal_Vein;
    Veinsignals2{find(VeinSliceNumbers==n_Veinslice)} = Signal_Vein2;
    Veinconcts{find(VeinSliceNumbers==n_Veinslice)} = Conct_Vein;
    Veincutoffs{find(VeinSliceNumbers==n_Veinslice)} = cutoffs_Vein;
    VeinsmoothFits{find(VeinSliceNumbers==n_Veinslice)} = SmoothVeinfit;
    VeinfitParams{find(VeinSliceNumbers==n_Veinslice)} = FitParam;
    Resnorms(find(VeinSliceNumbers==n_Veinslice)) = RESNORM;
    ibats(find(VeinSliceNumbers==n_Veinslice)) = ibat;
    bats(find(VeinSliceNumbers==n_Veinslice)) = bat;
    brts(find(VeinSliceNumbers==n_Veinslice)) = brt;
    tmptimes{find(VeinSliceNumbers==n_Veinslice)} = time;
tmptimes2{find(VeinSliceNumbers==n_Veinslice)} = time2;
tmpVeins{find(VeinSliceNumbers==n_Veinslice)} = Vein;
tmpVeins2{find(VeinSliceNumbers==n_Veinslice)} = Vein2;
tmpsmoothtimes{find(VeinSliceNumbers==n_Veinslice)} = time_smooth;
tmpsmoothVeins{find(VeinSliceNumbers==n_Veinslice)} = SmoothVeinfit;
        
%     if RESNORM < 1e-5
%         break
%     end

end
    
    %Loop to determine n_Veinslice
Resnorms_temp = Resnorms;
LoopingIndex = 1;
while LoopingIndex
    n_Veinslice = VeinSliceNumbers(find(Resnorms_temp==min(Resnorms_temp)));
    SmoothVeinfit = VeinsmoothFits{find(VeinSliceNumbers==n_Veinslice)};
    %Check to make sure n_Veinslice gives an actual Vein signal in shape
    if ( ( (max(SmoothVeinfit)-SmoothVeinfit(1)) > 0.015 ) && ( (SmoothVeinfit(length(SmoothVeinfit))-SmoothVeinfit(1)) < 0.015 ) )
        break;
    end  
    Resnorms_temp(find(VeinSliceNumbers==n_Veinslice)) = 999999;
    %Exit point
    if (min(Resnorms_temp) == 999999)
        n_Veinslice = n_Veinslice_max;
        break;
    end
end
    
    Mask_Vein_slice = VeinsliceMasks{find(VeinSliceNumbers==n_Veinslice)};
positionVein = Veinpositions{find(VeinSliceNumbers==n_Veinslice)};
Signal_Vein = Veinsignals{find(VeinSliceNumbers==n_Veinslice)};
Signal_Vein2 = Veinsignals2{find(VeinSliceNumbers==n_Veinslice)};
Conct_Vein = Veinconcts{find(VeinSliceNumbers==n_Veinslice)};
cutoffs_Vein = Veincutoffs{find(VeinSliceNumbers==n_Veinslice)};

FitParam = VeinfitParams{find(VeinSliceNumbers==n_Veinslice)};
RESNORM = Resnorms(find(VeinSliceNumbers==n_Veinslice));
ibat = ibats(find(VeinSliceNumbers==n_Veinslice));
bat = bats(find(VeinSliceNumbers==n_Veinslice));
brt = brts(find(VeinSliceNumbers==n_Veinslice));
SmoothVeinfit = VeinsmoothFits{find(VeinSliceNumbers==n_Veinslice)};

tmptime = tmptimes{find(VeinSliceNumbers==n_Veinslice)};
tmptime2 = tmptimes2{find(VeinSliceNumbers==n_Veinslice)};
tmpVein = tmpVeins{find(VeinSliceNumbers==n_Veinslice)};
tmpVein2 = tmpVeins2{find(VeinSliceNumbers==n_Veinslice)};
tmpsmoothtime = tmpsmoothtimes{find(VeinSliceNumbers==n_Veinslice)};
tmpsmoothVein = tmpsmoothVeins{find(VeinSliceNumbers==n_Veinslice)};

figure;plot(Signal_Vein);
title(['vein: ' num2str(bat) '-' num2str(brt) ' TR: ' num2str(TR/1000)]);
%exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_vein(auto).jpg'], 'resolution', 300);

figure();
plot(tmptime,tmpVein,'r*',tmptime2,tmpVein2,'k*',tmpsmoothtime,tmpsmoothVein,'k-');
xlabel('Time - sec');
ylabel('Gd Concentration');
title(['Vein - Smooth Fit - chosen slice ' num2str(n_Veinslice)]);
legend('Measured Vein','Gamma variate fit');
%exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_veinfit(auto).jpg'], 'resolution', 300);
%pause(0.1);

area = sum(SmoothVeinfit);
Vein = struct('Ct',SmoothVeinfit,'area',area,'FitParam',FitParam,'RESNORM',RESNORM,'ITP',ibat,'BATP',bat,'RTP',brt,'signal',Signal_Vein,'conc',Conct_Vein);
% Vein = struct('ITP',ibat,'BATP',bat,'Ct',Cmt,'RTP',brt,'area',area);
%     
% figure;
% subplot(2,1,1);
% plot(AIFslice.signal);title(['aif: ' num2str(AIFslice.BATP) '-' num2str(AIFslice.RTP)]);
% subplot(2,1,2);
% plot(Signal_Vein);hold on;plot(Signal_Vein2);
% title(['vein: ' num2str(bat) '-' num2str(brt) ' s' num2str(positionVein{1}(3)) ' id' caseid], 'interpreter', 'none');
% exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_aifautovein.jpg'], 'resolution', 300);

% subplot(2,2,3);
% plot(Mean_Signal);
% subplot(2,2,4);
% histogram(batlist,'binmethod','integers');title(caseid);

%yv(squeeze(dsc(:,:,1,positionVein{1}(3))), [], 'ol', Mask_Vein_slice(:,:,positionVein{1}(3)))
figure;imshow(squeeze(dsc(:,:,1,positionVein{1}(3))),[],'initialmagnification','fit');
    roir = ones(size(dsc,[1 2]));
    roib = zeros(size(dsc,[1 2]));
    roig = zeros(size(dsc,[1 2]));
    roic = cat(3, roir, roig, roib);
    hold on; imh = imshow(roic);
    set(imh, 'alphadata', 0.5*Mask_Vein_slice(:,:,positionVein{1}(3)));
    
%exportgraphics(gca, ['..\images\autovein debug\' currtime '_' caseid '_veinloc(auto).jpg'], 'resolution', 300);

end

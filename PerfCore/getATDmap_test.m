function [ATDmap,IBATmap,BATmap,BRTmap,sigmap] = getATDmap_test(AIF,dscstacks,automasks,path_DSC,Nmeas,Nslices,CmtWM,t1slc,Dt,VEIN)
%getATDmap Returns Arterial Tissue Delay map and voxel BAT map
%   all values are in data point number (not in sec or msec)
%   
% Author: Yong Ik Jeong
% Date: 2018-02-06

%% Get brain section ROIs
%roidata = load('..\SAVES\roiset_20180125v2.mat');
global glblTargetPath;
global injectionNum;
global seqType;

%     ATDmap = zeros(size(roistack));
%     BATmap = zeros(size(roistack));
%     IBATmap = zeros(size(roistack));
%     BRTmap = zeros(size(roistack));
%     sigmap = cell([length(unique(roistack))-1 2]);


% isvein = 0;
% if exist([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat'],'file')
%     % Get vein mask and vein signal
%     veindata = load([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat']);
%     veinmask = veindata.veinmask;
%     %path_LLEPI = [glblTargetPath '\P' sprintf('%03d',injectionNum) '\IR_LL_EPI_PRE'];
%     %path_DSC = [glblTargetPath '\P' sprintf('%03d',injectionNum) '\ep2d_perf'];
%     if size(veinmask,3) ~= 1 && size(veinmask,3) ~= Nslices
%         error('vein mask and dsc Nslices do not match');
%     end
%     if size(veinmask,3) ~= 1
%         [tmpi,tmpj,tmpz] = ind2sub(size(veinmask),find(veinmask,1));
%         veinmask = veinmask(:,:,tmpz);
%         if tmpz ~= t1slc
%             error('vein slc and t1 slc do not match');
%         end
%     end
%     vein_signal = zeros(1,Nmeas);
%     vein_signal2 = zeros(1,Nmeas);
%     for zz = 1:Nslices%size(veinmask,3)
%         dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(zz-1)*Nmeas),(zz*Nmeas),1)));
%         for ii = 1:size(veinmask,1)
%             for jj = 1:size(veinmask,2)
%                 if zz == t1slc && veinmask(ii,jj)
%                     vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),5,'sgolay',3),[1 Nmeas]);
%                     vein_signal2 = vein_signal2 + reshape(smooth(dscstack(ii,jj,:),.05,'sgolay',1),[1 Nmeas]);
%                     %vein_signal = vein_signal + reshape(dscstack(ii,jj,:),[1 Nmeas]);
%                     %vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),10),[1 Nmeas]);
%                 end
%             end
%         end
%     end
%     vein_signal = vein_signal / sum(veinmask(:));
%     vein_signal2 = vein_signal2 / sum(veinmask(:));
%     [vIBAT, vBAT, vBRT] = findBAT(vein_signal2);
%     [tmpibat,tmpbat,tmpbrt] = findBAT(vein_signal);
%     if abs(vBAT-tmpbat)*Dt > 1
%         if tmpbat > vBAT
%             vIBAT = tmpibat;
%             vBAT = tmpbat;
%             vBRT = tmpbrt;
%         end
%             
%     end
%         
%     isvein = 1;
%     figure;plot(vein_signal);hold on;plot(vein_signal2);
%     title(['vein: ' num2str(vBAT) '-' num2str(vBRT)]);
% end

isvein = 0;
if ~isempty(VEIN)
    isvein = 1;
    vBAT = VEIN.BATP;
    vBRT = VEIN.RTP;
end


% YIJ 20190507: get global average BAT
allavgsignal = zeros(1,Nmeas);
allsigcount = 0;
batlist = [];
brtlist = [];

%tmpfitopt = fitoptions('method','smoothingspline','smoothingparam',.1);
%tmpfittype = fittype('smoothingspline');
dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(1-1)*Nmeas),(1*Nmeas),1)));
smoothsignals = zeros(size(dscstack,1),size(dscstack,2),Nslices,Nmeas);
for zz = 1:Nslices
    if isempty(dscstacks)
        error('dscstacks is empty');
        dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(zz-1)*Nmeas),(zz*Nmeas),1)));
        tmpmask = automaskns(imfilter(dscstack,fspecial('gaussian',[3 3],2)));
        dscstack = fSpatialFilter(dscstack,[2 3],tmpmask);
    else
        dscstack = dscstacks{zz};
        tmpmask = automasks{zz};
    end
    
    for ii = 1:size(tmpmask,1)
        for jj = 1:size(tmpmask,2)
            if tmpmask(ii,jj)
                %tmpfit = fit([1:Nmeas]',squeeze(dscstack(ii,jj,:)),'smoothingspline','smoothingparam',.1);
                %tmpfit = fit([1:Nmeas]',squeeze(dscstack(ii,jj,:)),tmpfittype,tmpfitopt);
                %smoothsignal = reshape(tmpfit(1:Nmeas),[1 Nmeas]);
                %smoothsignal = reshape(smooth(dscstack(ii,jj,:),.05,'sgolay',1),[1 Nmeas]);
                
                % 20210326 YIJ: use 3rd order sgolay
                %smoothsignal = reshape(smooth(dscstack(ii,jj,:),max(sqrt(1/(Dt))*5,5),'sgolay',3),[1 Nmeas]);
                
                % 20210422 YIJ: remove smoothing because of arrival time shift
                smoothsignal = reshape(dscstack(ii,jj,:),[1 Nmeas]);
                % 20210329 YIJ: do spline smoothing
%                 tmpfit = fit([1:Nmeas]',squeeze(dscstack(ii,jj,:)),'smoothingspline','smoothingparam',.8);
%                 smoothsignal = reshape(tmpfit(1:Nmeas),[1 Nmeas]);
%                 smoothsignal(smoothsignal < 0) = 0;
                
                smoothsignals(ii,jj,zz,:) = smoothsignal;
                
                [tmpibat,tmpbat,tmpbrt] = findBAT_test(smoothsignal);
                if tmpbat > 5 && tmpbrt < 200
                    batlist = [batlist; tmpbat];
                    brtlist = [brtlist; tmpbrt];
                end
                tempsignal = smoothsignal;
                %tempsignal = reshape(dscstack(ii,jj,:),[1 Nmeas]);
                allavgsignal = allavgsignal + tempsignal;
                allsigcount = allsigcount + 1;
            end
        end
    end
end
allavgsignal = allavgsignal./allsigcount;
[allibat,allbat,allbrt] = findBAT_test(allavgsignal);
if allbat < AIF.BATP
    %error('global bat is earlier than aif bat');
end

fprintf('Global\t ibat %d | bat %d | brt %d\n', allibat, allbat, allbrt);
fprintf('AIF\t ibat %d | bat %d | brt %d\n', AIF.ITP, AIF.BATP, AIF.RTP);

if isempty(CmtWM)
    CmtWM = [allibat allbat allbrt];
end
if isempty(VEIN)
    vBAT = allbat;
    vBRT = allbrt;
end


ATDmap = zeros([size(tmpmask) Nslices]);
BATmap = zeros([size(tmpmask) Nslices]);
IBATmap = zeros([size(tmpmask) Nslices]);
BRTmap = zeros([size(tmpmask) Nslices]);
sigmap = {};
count = 0;
for zz = 1:Nslices%size(roistack,3)
%     if isempty(dscstacks)
%         dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(zz-1)*Nmeas),(zz*Nmeas),1)));
%         tmpmask = automaskns(imfilter(dscstack,fspecial('gaussian',[3 3],2)));
%         dscstack = fSpatialFilter(dscstack,[2 3],tmpmask);
%     else
%         dscstack = dscstacks{zz};
%     end
    
%     tmpallsig = reshape(dscstack,[size(dscstack,1)*size(dscstack,2) size(dscstack,3)]);
%     tmpmaskary = reshape(tmpmask,[size(tmpmask,1)*size(tmpmask,2) 1]);
%     tmpmasksig = mean(tmpallsig(logical(tmpmaskary),:),1);
%     [tmpibat,tmpbat,tmpbrt] = findBAT(tmpmasksig);
%     tmpmask = dscstack(:,:,1) > mean(tmpmasksig(tmpibat:tmpbat));
%     tmpmask = imfill(imdilate(tmpmask,strel('sphere',1)),'holes');

    tmpmask = automasks{zz};

    tmpbatmap = zeros(size(tmpmask));
    tmpibatmap = zeros(size(tmpmask));
    tmpbrtmap = zeros(size(tmpmask));
    tmpatdmap = zeros(size(tmpmask));

    for ii = 1:size(tmpmask,1)
        for jj = 1:size(tmpmask,2)
            if tmpmask(ii,jj)
                %if ii == 91 && jj == 44 && zz == 24
                if ii == 51 && jj == 53 && zz == 20
                    1;
                end
                signal = zeros(1,Nmeas);
                signal2 = zeros(1,Nmeas);
                
                tmpcount = 0;
                
                %tmpfit = fit([1:Nmeas]',squeeze(dscstack(ii,jj,:)),'smoothingspline','smoothingparam',.1);
%                 tmpfit = fit([1:Nmeas]',squeeze(dscstack(ii,jj,:)),tmpfittype,tmpfitopt);
%                 smoothsignal = reshape(tmpfit(1:Nmeas),[1 Nmeas]);
%                 %smoothsignal = reshape(smooth(dscstack(ii,jj,:),.05,'sgolay',1),[1 Nmeas]);
%                 smoothsignal(smoothsignal < 0) = 0;

                smoothsignal = reshape( smoothsignals(ii,jj,zz,:), [1 Nmeas]);

%                         if rr == ii && cc == jj
%                             signal = signal + weight(2)*reshape(smooth(dscstack(rr,cc,:),.05,'sgolay',1),[1 Nmeas]);
%                         else
%                             signal = signal + weight(1)*reshape(smooth(dscstack(rr,cc,:),.05,'sgolay',1),[1 Nmeas]);
%                         end
                        %signal = signal + reshape(smooth(dscstack(rr,cc,:),5,'sgolay',3),[1 Nmeas]);
                        %signal = signal + reshape(smooth(dscstack(rr,cc,:),.05,'sgolay',1),[1 Nmeas]);
                        %signal2 = signal2 + reshape(smooth(dscstack(rr,cc,:),5,'sgolay',1),[1 Nmeas]);
                        %signal2 = signal2 + reshape(dscstack(rr,cc,:),[1 Nmeas]);
                        %tmpcount = tmpcount + 1;

                %signal = signal ./ tmpcount;
                %signal2 = signal2 ./ tmpcount;
                %figure;plot(signal)
                [ibat,bat,brt] = findBAT_test(smoothsignal);
                % NEED TO ADD ALGORITHM FOR MULTIPLE PEAK (ARTERIAL AND
                % VENOUS) DETECTION
                if brt == 2
                    1;
                end
               
                % 20180305 YIJ: if bat and brt are greater than vein, set
                % to average wm bat and brt...
                if isvein
                    if bat > vBAT % 20180307 YIJ: +5 experimental
                        %figure;plot(signal);hold on;plot(signal2);
                        %title(['ROI: ' num2str(0) ' bat ' num2str(bat) '-' num2str(vBAT)]);
                        %bat = vBAT;
                        %bat = CmtWM(2);
                    end
                    if brt > vBRT
                        %figure;plot(signal);hold on;plot(signal2);
                        %title(['ROI: ' num2str(0) ' brt ' num2str(brt) '-' num2str(vBRT)]);
                        %brt = vBRT;
                        %brt = CmtWM(3);
                    end
                end
                % 20190411 YIJ: arrival time earlier than AIF by more than 1 sec is
                % probably wrong
                if (AIF.BATP - bat)*Dt > 1 
                    %bat = AIF.BATP;
                    %brt = AIF.RTP;
                    %ibat = AIF.ITP;
                    bat = AIF.BATP;
                    % If arrival time is likely wrong, set brt to AIF RTP instead)  
                    if (AIF.RTP - brt)*Dt > 1
                        %brt = AIF.RTP;
                        brt = allbrt;
                    end
                %end
                % 20190411 YIJ: arrival time later than vein BRT by more than 3 sec is
                % probably wrong
                elseif (bat - vBRT)*Dt > 5 %bat > vBRT
                %if (bat - vBAT)*Dt > 3
                    %bat = vBAT;
                    % If arrival time is likely wrong, set brt to vein BRT instead)
                    %brt = vBRT;
                %end
                
                % 20190417 YIJ: bolus width less than AIF width by more than 3 sec is
                % probably wrong so set to WM arrival times
                %if bat > AIF.BATP && bat < vBRT && (brt-bat)*Dt < (AIF.RTP-AIF.BATP)*Dt - 1
                elseif (brt-bat)*Dt < (AIF.RTP-AIF.BATP)*Dt - 4 &&...
                        (brt-bat)*Dt < mean(brtlist-batlist)*Dt - 4
                    %(brt-bat)*Dt < (AIF.RTP-AIF.BATP)*Dt - 3
                    %bat = CmtWM(2);
                    %brt = CmtWM(3);
                    %bat = AIF.BATP;
                    brt = allbrt;
                    
                %end
                
                % 20190424 YIJ: bolus width greater than vein and wm avg by
                % more than 10 sec is probably wrong so adjust to global
                % avg width
                %if bat > AIF.BATP && bat < vBRT && ((brt-bat)*Dt > (vBRT-vBAT)*Dt + 5 && (brt-bat)*Dt > (CmtWM(3)-CmtWM(2))*Dt + 5)
                elseif (brt-bat)*Dt > (vBRT-vBAT)*Dt + 10 && (brt-bat)*Dt > (CmtWM(3)-CmtWM(2))*Dt + 10 &&...
                        (brt-bat)*Dt > (allbrt-allbat)*Dt + 10
%                     bat = CmtWM(2);
%                     brt = CmtWM(3);
                    %brt = bat + (allbrt-allbat);
                    %brt = bat + max([allbrt-allbat CmtWM(3)-CmtWM(2) vBRT-vBAT]);
                    %brt = bat + round( mean(brtlist-batlist) + std(brtlist-batlist) );
                    brt = allbrt;
                    
                % 20210413 YIJ: If std dev changes after and before bolus
                % arrival are similar, then probably too noisy baseline
                % and/or too small of a bolus to detect, so set to AIF bat
                % and global brt
                elseif std(smoothsignal(5:brt))/std(smoothsignal(5:bat)) < 1.5
                    bat = AIF.BATP;
                    brt = allbrt;
                end
                
                if brt > length(smoothsignal)
                    brt = length(smoothsignal);
                end
                count = count + 1;
                tmpbatmap(ii,jj) = bat;
                tmpibatmap(ii,jj) = ibat;
                tmpbrtmap(ii,jj) = brt;
                tmpatdmap(ii,jj) = bat - AIF.BATP;
                sigmap{count,1} = count;
                sigmap{count,2} = [ii,jj,zz];
                %sigmap{count,3} = signal;
            end
        end
    end
    ATDmap(:,:,zz) = tmpatdmap;
    BATmap(:,:,zz) = tmpbatmap;
    IBATmap(:,:,zz) = tmpibatmap;
    BRTmap(:,:,zz) = tmpbrtmap;
end


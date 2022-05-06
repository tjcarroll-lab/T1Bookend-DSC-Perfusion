function [VEIN,veinmask,positionVein,n_Veinslice] = getVEIN(Dt,echoTime,path_DSC,Nmeas,Nslices,t1slc,n_slice_AIF,AIFslice)
%DDcorr Obtain Delay and Dispersion function fit parameters
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2019-06-27

global glblTargetPath;
global injectionNum;
global seqType;
global MANUAL;
global caseid;
global currtime;
if exist([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat'],'file')
    
%     caseid = split(path_DSC,{'\','/'});
%     caseid = caseid{end-3};
%     currtime = datestr(now,'yyyymmddHHMMSS');
    
    % Get vein mask and vein signal
    veindata = load([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat']);
    veinmask = veindata.veinmask;
    %path_LLEPI = [glblTargetPath '\P' sprintf('%03d',injectionNum) '\IR_LL_EPI_PRE'];
    %path_DSC = [glblTargetPath '\P' sprintf('%03d',injectionNum) '\ep2d_perf'];
    if size(veinmask,3) ~= 1 && size(veinmask,3) ~= Nslices
        error('vein mask and dsc Nslices do not match');
    end
    if size(veinmask,3) ~= 1
        [tmpi,tmpj,tmpz] = ind2sub(size(veinmask),find(veinmask,1));
        %veinmask = veinmask(:,:,tmpz);
        if isempty(t1slc)
            t1slc = tmpz;
        else
            if tmpz ~= t1slc
                error('vein slc and t1 slc do not match');
            end
        end
    end
    
    vein_signal = zeros(1,Nmeas);
    vein_signal2 = zeros(1,Nmeas);

    for s = 1:Nslices
%     tempimage = double(Cells2Matrix(IA_ReadImages(path_DSC,s,0,N_Slices)));
        tempimage = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(s-1)*Nmeas),(s*Nmeas),1))); 
        Images{s} = tempimage;
    end

    % 20210326 YIJ: spatially smooth DSC (xy and z)
    Images = smoothDSC(Images);

    %dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(t1slc-1)*Nmeas),(t1slc*Nmeas),1)));
    dscstack = Images{t1slc};
    for ii = 1:size(veinmask,1)
        for jj = 1:size(veinmask,2)
            if veinmask(ii,jj,t1slc)
                %vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),5,'sgolay',3),[1 Nmeas]);
                %vein_signal2 = vein_signal2 + reshape(smooth(dscstack(ii,jj,:),max(sqrt(1/(Dt))*5,5),'sgolay',3),[1 Nmeas]);
                %vein_signal = vein_signal + reshape(dscstack(ii,jj,:),[1 Nmeas]);
                %vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),10),[1 Nmeas]);
                
                % 20210422 YIJ: remove smoothing because of arrival time shift
                vein_signal = vein_signal + reshape(dscstack(ii,jj,:),[1 Nmeas]);
                % 20210329 YIJ: do spline smoothing
%                 tmpfit = fit([1:Nmeas]',squeeze(dscstack(ii,jj,:)),'smoothingspline','smoothingparam',.8);
%                 vein_signal = vein_signal + reshape(tmpfit(1:Nmeas),[1 Nmeas]);
%                 vein_signal(vein_signal < 0) = 0;
                vein_signal2 = vein_signal;

            end
        end
    end
    
    vein_signal = vein_signal / sum(veinmask(:));
    vein_signal2 = vein_signal2 / sum(veinmask(:));
    
%     vein_signal = zeros(1,Nmeas);
%     vein_signal2 = zeros(1,Nmeas);
%     for zz = 1:Nslices%size(veinmask,3)
%         dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(zz-1)*Nmeas),(zz*Nmeas),1)));
%         for ii = 1:size(veinmask,1)
%             for jj = 1:size(veinmask,2)
%                 if zz == t1slc && veinmask(ii,jj)
%                     vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),5,'sgolay',3),[1 Nmeas]);
%                     vein_signal2 = vein_signal2 + reshape(smooth(dscstack(ii,jj,:),.05,'sgolay',3),[1 Nmeas]);
%                     %vein_signal = vein_signal + reshape(dscstack(ii,jj,:),[1 Nmeas]);
%                     %vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),10),[1 Nmeas]);
%                 end
%             end
%         end
%     end
%     vein_signal = vein_signal / sum(veinmask(:));
%     vein_signal2 = vein_signal2 / sum(veinmask(:));
    
    [vIBAT, vBAT, vBRT] = findBAT_test(vein_signal2);
    [tmpibat,tmpbat,tmpbrt] = findBAT_test(vein_signal);
    if abs(vBAT-tmpbat)*Dt > 1
        if tmpbat > vBAT
            vIBAT = tmpibat;
            vBAT = tmpbat;
            vBRT = tmpbrt;
        end
            
    end
    vbats = [vIBAT vBAT vBRT];
    VEIN = Sig2Conct_AIF(vein_signal,vbats,echoTime);
    VEIN.signal = vein_signal;
    VEIN.conc = VEIN.Ct;
    VEIN.Ct = VEIN.Ct(VEIN.BATP:VEIN.RTP);
    
    figure;plot(vein_signal);hold on;plot(vein_signal2);
    title(['vein: ' num2str(vBAT) '-' num2str(vBRT) ' TR: ' num2str(Dt)]);
    %exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_vein(manu).jpg'], 'resolution', 300);
    
    
    % vein t1 enhancement correction
    %[VEIN.Ct,VEIN.RTP] = enhanceCorr(AIF,vein_signal,echoTime);
   
    % fit vein signal to gamma variate
    vtime = 0:Dt:(VEIN.RTP-VEIN.BATP)*Dt;
    BestFit = gamma_variate_IDC_gridsearch(vtime,VEIN.Ct);
    SumSquares = BestFit(1);
    A = BestFit(2);
    alpha = BestFit(3);
    beta = BestFit(4);
    [vA,valpha,vbeta,RESNORM,RESIDUAL,EXITFLAG] = gamma_variate_IDC_LeastSquaresFit(vtime,VEIN.Ct,A,alpha,beta);
    FitParam.A = vA;
    FitParam.alpha = valpha;
    FitParam.beta = vbeta;
    time_smooth = 0:Dt/Nslices:(VEIN.RTP-VEIN.BATP)*Dt;
    tmpvein = gamma_variate_model_IDC(time_smooth,vA,valpha,vbeta);
    figure();
    plot(vtime,VEIN.Ct,'r*',time_smooth,tmpvein,'k-');
    xlabel('Time - sec');
    ylabel('Gd Concentration');
    title('Vein - Smooth Fit');
    legend('Measured Vein','Gamma variate fit');
    %exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_veinfit(manu).jpg'], 'resolution', 300);
    %figure;plot(vtime,VEIN.Ct(VEIN.BATP:VEIN.RTP));hold on;plot(vtime,tmpvein);
    VEIN.Ct = tmpvein;
    VEIN.FitParam = FitParam;
    
    count = 1;
    clear positionVein
    for i = 1:size(veinmask,1)
        for j = 1:size(veinmask,2)
            if veinmask(i,j,t1slc)
                positionVein{count} = [i,j,t1slc];
                count = count + 1;
            end
        end
    end
    n_Veinslice = t1slc;
    VEIN.auto = 0;
    
%     figure;
% subplot(2,1,1);
% plot(AIFslice.signal);title(['aif: ' num2str(AIFslice.BATP) '-' num2str(AIFslice.RTP)]);
% subplot(2,1,2);
% plot(vein_signal);hold on;plot(vein_signal2);
% title(['vein: ' num2str(vBAT) '-' num2str(vBRT) ' s' num2str(positionVein{1}(3)) ' id' caseid], 'interpreter', 'none');
% exportgraphics(gcf, ['..\images\autovein debug\' currtime '_' caseid '_aifmanuvein.jpg'], 'resolution', 300);

% subplot(2,2,3);
% plot(Mean_Signal);
% subplot(2,2,4);
% histogram(batlist,'binmethod','integers');title(caseid);

%yv(dscstack(:,:,1), [], 'ol', veinmask)
figure;imshow(dscstack(:,:,1),[],'initialmagnification','fit');
    roir = ones(size(dscstack,[1 2]));
    roib = zeros(size(dscstack,[1 2]));
    roig = zeros(size(dscstack,[1 2]));
    roic = cat(3, roir, roig, roib);
    hold on; imh = imshow(roic);
    set(imh, 'alphadata', 0.5*veinmask(:,:,t1slc));
    
%exportgraphics(gca, ['..\images\autovein debug\' currtime '_' caseid '_veinloc(manu).jpg'], 'resolution', 300);

else
    if MANUAL
        error('Cant find VEIN mask for finding vein');
    else
        warning('Cant find VEIN mask for finding vein; using auto instead');
        [VEIN,veinmask,positionVein,n_Veinslice] = autoVein(path_DSC,Nmeas,Nslices,AIFslice);
        VEIN.auto = 1;
    end
end
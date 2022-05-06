function [DDfit,VEIN] = DDcorr_test(AIF,Dt,echoTime,Wij,path_DSC,Nmeas,Nslices,t1slc)
%DDcorr Obtain Delay and Dispersion function fit parameters
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2018-02-05

%% AIF and VEIN deconvolution

global glblTargetPath;
global injectionNum;
global seqType;
if exist([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat'],'file')
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
        veinmask = veinmask(:,:,tmpz);
        if tmpz ~= t1slc
            error('vein slc and t1 slc do not match');
        end
    end
    
    vein_signal = zeros(1,Nmeas);
    vein_signal2 = zeros(1,Nmeas);
    for zz = 1:Nslices%size(veinmask,3)
        dscstack = double(Cells2Matrix(IA_ReadImages(path_DSC,(1+(zz-1)*Nmeas),(zz*Nmeas),1)));
        for ii = 1:size(veinmask,1)
            for jj = 1:size(veinmask,2)
                if zz == t1slc && veinmask(ii,jj)
                    vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),5,'sgolay',3),[1 Nmeas]);
                    vein_signal2 = vein_signal2 + reshape(smooth(dscstack(ii,jj,:),.05,'sgolay',1),[1 Nmeas]);
                    %vein_signal = vein_signal + reshape(dscstack(ii,jj,:),[1 Nmeas]);
                    %vein_signal = vein_signal + reshape(smooth(dscstack(ii,jj,:),10),[1 Nmeas]);
                end
            end
        end
    end
    vein_signal = vein_signal / sum(veinmask(:));
    vein_signal2 = vein_signal2 / sum(veinmask(:));
    
    [vIBAT, vBAT, vBRT] = findBAT(vein_signal2);
    [tmpibat,tmpbat,tmpbrt] = findBAT(vein_signal);
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
    tmpvein = gamma_variate_model_IDC(vtime,vA,valpha,vbeta);
    %figure;plot(vtime,VEIN.Ct(VEIN.BATP:VEIN.RTP));hold on;plot(vtime,tmpvein);
    VEIN.Ct = tmpvein;
    VEIN.FitParam = FitParam;
    
    atime = 0:Dt:(AIF.RTP-AIF.BATP)*Dt;
    aA = AIF.FitParam.A;
    aalpha = AIF.FitParam.alpha;
    abeta = AIF.FitParam.beta;
    
    % set a threshold to cut out front-end zeros (e.g. inaccurate BAT)
    [~,nmax] = max(AIF.Ct);
    cutoffind = find(AIF.Ct(1:nmax) > 0.0001,1);
    AIF.Ct = AIF.Ct(cutoffind - 1:end);
    AIF.BATP = AIF.BATP + cutoffind - 1 - 1;
    
    [~,nmax] = max(VEIN.Ct);
    cutoffind = find(VEIN.Ct(1:nmax) > 0.0001,1);
    VEIN.Ct = VEIN.Ct(cutoffind - 1:end);
    VEIN.BATP = VEIN.BATP + cutoffind - 1 - 1;
    
    % smoothing AIF
    N = length(AIF.Ct);
    input_AIF = [AIF.Ct(1) (AIF.Ct(1:N-2)+4*AIF.Ct(2:N-1)+AIF.Ct(3:N))/6 AIF.Ct(N)];
    windowdiff = length(VEIN.BATP:VEIN.RTP) - length(AIF.BATP:AIF.RTP);
    if windowdiff < 0
        addaif = 0;
        addct = abs(windowdiff);
    else
        addaif = windowdiff;
        addct = 0;
    end
    
    % add extra points to make AIF and VEIN equal length
    input_AIF = [input_AIF gamma_variate_model_IDC(atime(end)+Dt:Dt:atime(end) + addaif*Dt,aA,aalpha,abeta)];
    input_VEIN = [VEIN.Ct gamma_variate_model_IDC(vtime(end)+Dt:Dt:vtime(end) + addct*Dt,vA,valpha,vbeta)];
    
    n = length(input_AIF);
    
    % make AIF matrix
    for i = 1:n
        for j=1:i
            AIF_Mat(i,j) = Dt*input_AIF(i-j+1) ;
        end
    end
    
    % sigular value decomposition
    [U,S,V] = svd(AIF_Mat);
    inv_S = invert_S(S,Wij);
    invA_mat = V*inv_S*U';
    
    % calculation Residue function
    Residue = (invA_mat*input_VEIN')';
    %figure;plot(Residue);
    [~,nmax] = max(Residue);
    Residue = Residue(nmax:end);
    
    %% Venous residue fitting
    
    ATD = Dt*abs(VEIN.BATP - AIF.BATP);
    
    fiteq = fittype('a/(ATD + 1)*exp(-b*x/ATD)','problem','ATD');
    time = 0:Dt:(length(Residue)-1)*Dt;
    [fit0,gof0,output0] = fit(time',Residue',fiteq,'problem',ATD);
    DDfit = fit0;
    
    figure;
    subplot(1,2,1);plot(input_AIF);hold on;plot(input_VEIN);
    subplot(1,2,2);plot(time,Residue);hold on;plot(fit0);
    
    
    %%
    % ddfunc = (fit0.a./(ATD + 1).*exp(-fit0.b.*[0:Dt:(size(AIF_Mat,1)-1)*Dt]/ATD))';
    % newAIF = AIF_Mat*ddfunc;
    
else
    DDfit = [];
end
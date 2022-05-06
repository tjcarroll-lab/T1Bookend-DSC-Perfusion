function [rCBV_DSC,rCBF_SVD,CMTT_CVP,Tmax_map,TTP_map,dBATmap,Conc_map] = fCBV_CBF_MTT_nSVD_DD(images,mask,AIF,Dt,Wij, echoTime,CmtWM,DDfit,ATDmap,IBATmap,BATmap,BRTmap)
% Mira change above code to make sure Tmax_map and TTP_map are outputs-----
%%=======================================================================%%
%%          rCBV,rCBF,MTT,1/SNR_c calculation function         
%%=======================================================================%%
% input :   1) images (ixj matrix)  : DSC images
%           2) mask (ixj martrix)   : non-sculp maks
%           3) AIF (structure)
%               AIF.Ct (1xN vector) : Conct. time data
%               AIF.BATP (scalar)   : Bolus arrival time for AIF
%           4) Dt (scalar) : discrete time interval(TR)
%           5) Wij (scalar)         : Threshold in SVD
%%=======================================================================%%
% output:   1) CBF_SVD  (2-D martrix) : rCBF in SVD
%           2) CMMT_SVD (2-D martrix) : CMTT from Residue fuction in SVD
%           3) CBV_CVP  (2-D martrix) : = CBF_SVD * CMMT_SVD 
%           4) CBV_DSC  (2-D martrix) : rCBV from the ratio of integral of
%           5) CMTT_CVP (2-D martrix) : = CBV_SDSC / CBF_SVD 
%           6) reciSNR_c(2-D martrix) : 3-D SNR in Conct.
%                                                               02-262-2004 

% define variables
Kh = (1-0.45)/(1-0.25);
psi = 1.04; % mg/ml --> g/mL?
global sampdelayfix;
global DDfix;

% using 10 points to calculate CBV
% if (AIF.BATP-AIF.RTP) <=10
%     if (AIF.BATP+9)  <= length(AIF.Ct)
%         AIF.area = sum(AIF.Ct(AIF.BATP:AIF.BATP+9));
%     else
%         AIF.area = sum(AIF.Ct(AIF.BATP:length(AIF.Ct)));
%     end
% else
%     AIF.area = sum(AIF.Ct(AIF.BATP:AIF.RTP));
% end

% initial value in map
rCBV_DSC  = zeros(size(images,1), size(images,2));
rCBF_SVD  = zeros(size(images,1), size(images,2));
CMTT_CVP  = zeros(size(images,1), size(images,2));
% Mira change this code----------------------------------------------------
Tmax_map   = zeros(size(images,1), size(images,2));
dBATmap   = zeros(size(images,1), size(images,2));
Conc_map  = zeros(size(images,1), size(images,2));
TTP_map = zeros(size(images,1), size(images,2));
% -------------------------------------------------------------------------

if isfield(AIF,'ITP')
    ITP=AIF.ITP; BAT=AIF.BATP; RTP=AIF.RTP;
else
    ITP=2; BAT=AIF.BATP; RTP=AIF.RTP;
end

warning off;

% save original AIF
oAIF = AIF;

%YIJ 20170616 parallel computing
jsize = size(images,2);
% if ~matlabpool('size')
%     matlabpool
% end
%parfor i = 1:size(images,1)

tmpfitopt = fitoptions('method','smoothingspline','smoothingparam',.1);
tmpfittype = fittype('smoothingspline');
tmpcount = 0;
tmpdebug = 0;

for i = 1:size(images,1)
    for j = 1:jsize
        if mask(i,j) && BATmap(i,j) %YIJ 20190403: BATmap has background masked out (only for getATDmap_test.m)
            
            if BRTmap(i,j) ~= 0 %%ATDmap(i,j) == 7 %&& BRTmap(i,j) == 70
                1;
            end
                      
            % call info. from Cmt strucure
            Signal = squeeze(images(i,j,:))';
            
            if (i == 100 && j == 100) %|| (i == 133 && j == 92) %|| (i == 55 && j == 108) %20190718 YIJ: test for case 23
                1;
            end
            
            oSignal = Signal;
            
            %YIJ 20170613         
            tmpfit = fit([1:length(Signal)]',Signal',tmpfittype,tmpfitopt);
            Signal = reshape(tmpfit(1:length(Signal)),[1 length(Signal)]);
            Signal(Signal < 0) = 0;
            %Signal = smooth(Signal,5,'sgolay',3)';
            %Signal = smooth(Signal,.05,'sgolay',1)';
            %Signal = freqfilter(Signal,'low',10);
            %[ccv,vBET] = enhanceCorr(AIF,Signal,echoTime);
            
            [isbolus,scndfilt] = boluscheck(Signal,oSignal,oAIF,[IBATmap(i,j) BATmap(i,j) BRTmap(i,j)], echoTime);
            tmpdebug = tmpdebug + scndfilt;
            tmpcount = tmpcount + 1;
                                   
            if DDfix && BATmap(i,j) && IBATmap(i,j) && BRTmap(i,j)
                % get regional conc. curve
                Cmt = Sig2Conct_AIF(Signal, [IBATmap(i,j) BATmap(i,j) BRTmap(i,j)], echoTime);
                Cmt2 = Sig2Conct_AIF(Signal, [IBATmap(i,j) BAT BRTmap(i,j)], echoTime);
                % apply DDfix to AIF
                AIF = applyDDfix_test(oAIF,Dt,ATDmap(i,j),DDfit,Wij);
                %tmp = AIF;
                % Fix AIF ratio
                AIF.Ct = AIF.Ct .* sum(oAIF.Ct)/sum(AIF.Ct);
                if abs(sum(oAIF.Ct) - sum(AIF.Ct)) > .0000001
                    1;
                end
                if ATDmap(i,j) > 5
                    1;
                end
            else
                %Cmt = Sig2Conct_AIF(Signal, [CmtWM(1:3)], echoTime);
                Cmt = Sig2Conct_AIF(Signal, [IBATmap(i,j) BAT BRTmap(i,j)], echoTime);
                Cmt2 = Sig2Conct_AIF(Signal, [IBATmap(i,j) BAT length(Signal)], echoTime);
%                 Signal = smooth(Signal,.05,'sgolay',1)';
%                 Cmt2 = Sig2Conct_AIF(Signal, [IBATmap(i,j) BAT BRTmap(i,j)], echoTime);
                AIF = oAIF;
            end
            if Cmt.RTP == 0
                1;
            end
            % calculate Residue function
            %fResidue = fResidue_nSVD(AIF,Cmt.Ct,Dt,Wij);
            fResidue = fResidue_nSVD(AIF,Cmt,Dt,Wij);
            
            fResidue2 = fResidue_nSVD(AIF,Cmt2,Dt,Wij);
            
%             if BATmap(i,j) && max(fResidue) < 0.01
%                 disp(i);
%                 disp(j);
%             end
            % make Tmax map
            [MaxIn MaxN] = max(fResidue);
            [MaxIn2 MaxN2] = max(fResidue2);
            
            if max(fResidue) == 0
                1;
            end
            
            % store info
            %reciSNR_c(i,j) = Cmt.reciSNR_c;
            %dBATmap(i,j)   = Dt*(Cmt.BATP-AIF.BATP); %NOTE (Jessy - 01/05/10): This is using same BATP for tissue as AIF (needs new computation of BATP of tissue curve using "findBAT" function)
            %dBATmap(i,j)   = Dt*(Cmt.BATP-AIF.BATP);
            dBATmap(i,j) = isbolus;
            rCBF_SVD(i,j)  = 100*60*Kh/psi*max(fResidue);  % ml/100mg/min 
            %rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(AIF.BATP:AIF.RTP))/sum(AIF.Ct(AIF.BATP:AIF.RTP)); % ml/100g
            %rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(Cmt.BATP:Cmt.RTP))/sum(AIF.Ct(AIF.BATP:AIF.RTP));
            if sampdelayfix
                % YIJ 201904 Try using original (no DD fix) AIF to keep CBV
                % same
                rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(Cmt.BATP:Cmt.RTP))/sum(AIF.Ct);
            else
                rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(Cmt.BATP:Cmt.RTP))/sum(AIF.Ct(AIF.BATP:AIF.RTP));
            end
            % Mira Add this code-------------------------------------------
            Tmax_map(i,j)   = Dt*(MaxN-1);
            [~,peakind] = max(Cmt.Ct);
            TTP_map(i,j) = Dt*(peakind-1);
            % -------------------------------------------------------------
            
            Conc_map(i,j)  = Dt*(MaxN2-1);
%             if isbolus
%                 TmaxMap(i,j)   = Dt*(MaxN-1);
%                 %Conc_map(i,j)  = Dt*(MaxN2-1);
%             else
%                 TmaxMap(i,j) = Inf;
%                 %Conc_map(i,j)  = Inf;
%             end           
%             Conc_map(i,j)  = Dt*(MaxN-1);
            CMTT_CVP(i,j)  = 60*rCBV_DSC(i,j)/rCBF_SVD(i,j); %sec
            %Conc_map(i,j)  = sum(Cmt.Ct(Cmt.BATP:Cmt.RTP));
            %Conc_map(i,j)  = sum(AIF.Ct)/sum(oAIF.Ct);
            
        end
    end
end

% YIJ 20200113: bolus filtering debug
fprintf('bolus secondary filter rate: %d / %d = %f\n',tmpdebug,tmpcount,tmpdebug/tmpcount);

% clear variables
clear images mask
function [rCBV_DSC,rCBF_SVD,CMTT_CVP,TmaxMap,dBATmap,Conc_map] = fCBV_CBF_MTT_nSVD(images,mask,AIF,Dt,Wij, echoTime,CmtWM)

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
psi = 1.04; % mg/ml
global sampdelayfix;

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
TmaxMap   = zeros(size(images,1), size(images,2));
dBATmap   = zeros(size(images,1), size(images,2));
Conc_map  = zeros(size(images,1), size(images,2));

if isfield(AIF,'ITP')
    ITP=AIF.ITP; BAT=AIF.BATP; RTP=AIF.RTP;
else
    ITP=2; BAT=AIF.BATP; RTP=AIF.RTP;
end

warning off;

%YIJ 20170616 parallel computing
jsize = size(images,2);
% if ~matlabpool('size')
%     matlabpool
% end
%parfor i = 1:size(images,1)
for i = 1:size(images,1)
    for j = 1:jsize
        if mask(i,j) 
            
            if (j == 115 && i == 60)
                1;
            end
            
            
            % call info. from Cmt strucure
            Signal = squeeze(images(i,j,:))';
            
            %YIJ 20170613
            Signal = smooth(Signal,5,'sgolay',3)';
            %Signal = freqfilter(Signal,'low',10);
            
            %Cmt = Sig2Conct_AIF(Signal,[ITP BAT 200], echoTime); %NOTE (Jessy - 01/05/10): should give individually computed BAT for tissue as an argument)           
            Cmt = Sig2Conct_AIF(Signal, CmtWM(1:3), echoTime);
            %Cmt = Sig2Conct_AIF(Signal,[ITP BAT CmtWM(3)], echoTime);
            %Cmt = Sig2Conct_AIF(Signal,[CmtWM(1) BAT CmtWM(3)], echoTime);
            
            % calculate Residue function
            %fResidue = fResidue_nSVD(AIF,Cmt.Ct,Dt,Wij);
            fResidue = fResidue_nSVD(AIF,Cmt,Dt,Wij);
            
            % make Tmax map
            [MaxIn MaxN] = max(fResidue);
            
            % store info
            %reciSNR_c(i,j) = Cmt.reciSNR_c;
            %dBATmap(i,j)   = Dt*(Cmt.BATP-AIF.BATP); %NOTE (Jessy - 01/05/10): This is using same BATP for tissue as AIF (needs new computation of BATP of tissue curve using "findBAT" function)
            dBATmap(i,j)   = Dt*(Cmt.BATP-AIF.BATP);
            rCBF_SVD(i,j)  = 100*60*Kh/psi*max(fResidue);  % ml/100mg/min 
            %rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(AIF.BATP:AIF.RTP))/sum(AIF.Ct(AIF.BATP:AIF.RTP)); % ml/100g
            %rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(Cmt.BATP:Cmt.RTP))/sum(AIF.Ct(AIF.BATP:AIF.RTP));
            if sampdelayfix
                rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(Cmt.BATP:Cmt.RTP))/sum(AIF.Ct);
            else
                rCBV_DSC(i,j)  = 100*Kh/psi*sum(Cmt.Ct(Cmt.BATP:Cmt.RTP))/sum(AIF.Ct(AIF.BATP:AIF.RTP));
            end
            TmaxMap(i,j)   = Dt*MaxN;
            CMTT_CVP(i,j)  = 60*rCBV_DSC(i,j)/rCBF_SVD(i,j); %sec
            Conc_map(i,j)  = sum(Cmt.Ct(Cmt.BATP:Cmt.RTP));
        end
    end
end

% clear variables
clear images mask
function AIF = Sig2Conct_AIF(Sig_AIF,cutoff,TE)

%%=======================================================================%%
%%          Convert signal into concentation for AIF        
%%=======================================================================%%
% input :   1) Sig_AIF (1xN vector) : AIF signal intensity data
%           2) cutoff ([i j])       : estimated BATP and RTP in AIF
%%=======================================================================%%
% output:   1) AIF (structure)
%               AIF.Ct (1xN vector) : Conct. time data
%               AIF.ITP (scalar)    : starting time point in AIF
%               AIF.BATP (scalar)   : Bolus arrival time point in AIF
%               AIF.RTP (scalar)    : Recirculation time point in AIF
%               AIF.PTP (scalar)    : Peak time point in AIF
%                                                               01-09-2007 


if nargin < 3
 %   TE = 1;
end

if length(cutoff) < 3
    sBAT = 2; iBAT=cutoff(1);   rBAT = cutoff(2);
else
    sBAT = cutoff(1); iBAT=cutoff(2);   rBAT = cutoff(3);
end

% calculate initial signal by approximate BAT ~ cutoff(1)
AverS0 = mean(Sig_AIF(sBAT:iBAT)); % Averaged before first cutoff time point except first one

% convert signal into concentration
Cmt = -log(Sig_AIF./AverS0)/TE;
Cmt(find(~isfinite(Cmt)))=0;

%Grady addition 10_26_2016
Cmt(find(Cmt < 0))=0;

% area
area = sum(Cmt(iBAT:rBAT));

% % find maximum signal intensity 
% [In_max,m] = max(temp_Cmt(iBAT:rBAT)); 
% BATP = cutoff(1);

% % find BATP to make more than 10% change to campare with maximum Intensity
% count = cutoff(1);
% for i = cutoff(1):length(Sig_AIF)
%     temp_ratio = temp_Cmt(i)/In_max;
%     if temp_ratio > 0.10
%         BATP = count;
%         break;
%     end;
%     count = count + 1;
% end;

% if temp_AverS0 
% % SNR calculation
% Smin = min(Sig_AIF);
% Sigma0 = std(Sig_AIF(2:BATP-1));
% sigma = sqrt(BATP)*Sigma0;
% AverS0 = mean(Sig_AIF(2:BATP-1));
% 
% % store info
% Ct   = -1/TE*log(Sig_AIF./AverS0); 
% Ct(find(~isfinite(Ct)))=0;
% RTP  = cutoff(2);
% area = sum(Ct(BATP:RTP));
% 
% else
%     % store info
% Ct   = zeros(1,length(Sig_AIF));
% BATP = cutoff(1); 
% RTP  = cutoff(2);
% area = 0;
% %AIF.reciSNR_c  = 0;
% end

AIF = struct('Ct',Cmt,'ITP',sBAT,'BATP',iBAT,'RTP',rBAT,'area',area,'AvgS0',AverS0);
% clear Sig_AIF cutoff temp_AverS0 temp_Cmt In_max m temp_ratio count Smin Sigma0 sigma AverS0 temp BATP Ct RTP area
function [AIF] = applyDDfix_test(AIF,Dt,ATD,DDfit,Wij)
%applyDDfix apply delay and dispersion fix to AIF
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2018-02-07

% apply DDfix to AIF
%if ATD > 0
% YIJ 20200826
if ATD*Dt > 1
    AIF.Ct = [AIF.Ct zeros(1,ATD)];
end
N = length(AIF.Ct);
%input_AIF = [AIF.Ct(1) (AIF.Ct(1:N-2)+4*AIF.Ct(2:N-1)+AIF.Ct(3:N))/6 AIF.Ct(N)];
input_AIF = AIF.Ct;
for r = 1:N
    for c=1:r
        AIF_Mat(r,c) = Dt*input_AIF(r-c+1);
    end
end
ATDt = abs(ATD)*Dt;
DDfun = (DDfit.a/(ATDt + 1)).*exp(-DDfit.b.*[0:Dt:(N-1)*Dt]/ATDt);
%DDfun = (DDfit.b/(ATDt)).*exp(-DDfit.b.*[0:Dt:(N-1)*Dt]/ATDt);
for r = 1:N
    for c=1:r
        DD_Mat(r,c) = Dt*DDfun(r-c+1);
    end
end

% skip ATD -3 to +3 since it could amplify noise
%if ATD > 0
% YIJ 20200826
if ATD*Dt > 1
    % convolve AIF with dispersion func
    %AIF.Ct = (AIF_Mat*DDfun')';
    AIF.Ct = conv(input_AIF,DDfun)*Dt;
    AIF.oRTP = AIF.RTP;
    %AIF.RTP = AIF.RTP + ATD;
    AIF.RTP = AIF.RTP + ATD + length(DDfun) - 1;
%elseif ATD < -5
% YIJ 20200826
elseif ATD*Dt < -1
    % deconvolve
    [U,S,V] = svd(DD_Mat);
    inv_S = invert_S(S,Wij);
    inv_mat = V*inv_S*U';
    AIF.Ct = (inv_mat*input_AIF')';
else
end


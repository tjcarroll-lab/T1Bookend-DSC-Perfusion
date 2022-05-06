function [ccv,vBET] = enhanceCorr(AIF,vein,echoTime)
%enhanceCorr correct for T1 enhancement in vein
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2018-02-23

[vIBAT vBAT vBRT] = findBAT(vein);
vbats = [vIBAT vBAT vBRT];
[a,b,c] = findBAT(vein(end:-1:1));
vBET = length(vein)-b+1;

% convert to conc.
VEIN = -1/echoTime*log(vein./mean(vein(vIBAT:vBAT)));
VEIN(find(~isfinite(VEIN)))=0;

ca = AIF.Ct;
cv = VEIN(vBAT:vBET);

%figure;plot(ca);hold on;plot(cv);

% leakage profile?
integral = zeros(size(ca));
for ii = 1:length(ca)
    integral(ii) = sum(ca(1:ii));
end

% leakage factor
lf = abs(mean(VEIN(vBET:end))) / max(integral); 

integral2 = interp1([1:length(ca)],integral,linspace(1,length(ca),length(cv)));

% corrected vein signal
ccv = cv + lf*integral2;
%figure;plot(cv);hold on;plot(ccv);

%figure;plot(ca);hold on;plot(cv);plot(ccv);plot(lf*integral2);
end


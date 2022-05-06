function [ROIsMaskWM,Threshold,T1distribution,T1table] = fAutoMaskWM_SS(varargin)
%(T1map,Tesla,mask,dT1)

% Get automask for white matter from T1 map
% function : ROIsMaskWM = fAutoMaskWM(T1map,Tesla,mask)
%
% Authoer : Wanyong Shin
%
% description: WM mask is acquired based on T1 distribution 
% using half-hight cutting off.
%             
% status: stable
%
% versions
%   02-09-2005 (WYS): Initialization
%   05-18-2005 (WYS): updated
%   08-09-2006 (JJM): updated

if nargin < 1;  sprint('more than one input variable is required');
else        ;   T1map = varargin{1};    end
if nargin < 2;  Tesla = '';
else        ;   Tesla = varargin{2};    end
if nargin < 3;  mask = AutoMaskNSforWM(T1map);
else         ;  mask = varargin{3}; end
if nargin < 4;  dT1 = 20;
else         ;  dT1 = varargin{4};  end
if nargin < 5;  cutNpixel = 400;
else         ;  cutNpixel = varargin{5};  end

if strmatch ('3.0T',Tesla)
    appoxHighT1WM = 750;
else
    appoxHighT1WM = 600;
end

% definition of variables
tempthreshold = dT1*round(200/dT1);
T1table = 0:dT1:2000;

% get T1 distribution from T1 map
T1distribution = hist(T1map(find(mask)),T1table);

% calculate Maximum T1 value in the range of appoxHighT1WM +/- tempthreshold
[temp,Index] = min(abs(T1table-appoxHighT1WM));
[nMaxT1WM,b] = max(T1distribution(Index-tempthreshold/dT1:Index+tempthreshold/dT1));
nIndex = Index-tempthreshold/dT1+b-1;   % 
MaximumT1WM = T1table(nIndex);      % maximum T1 value (600ms)

% find high/low threshold with half hight of maximum T1 from T1 distribution 
Lcount = nIndex;    Hcount = nIndex;
while T1distribution(Hcount) - nMaxT1WM/2 > 0
%     HighThreshold = T1table(count);
    Hcount = Hcount + 1;
end
while T1distribution(Lcount) - nMaxT1WM/2 > 0
%     LowThreshold = T1table(count);
    Lcount = Lcount - 1;
end

% calcualte threshold
y1 = T1distribution(Lcount);
y2 = T1distribution(Lcount+1);
LowThreshold = round(T1table(Lcount) + (nMaxT1WM/2-y1)/(y2-y1)*(T1table(Lcount+1) - T1table(Lcount)));

y1 = T1distribution(Hcount-1);
y2 = T1distribution(Hcount);
HighThreshold = round(T1table(Hcount-1) + (nMaxT1WM/2-y1)/(y2-y1)*(T1table(Hcount) - T1table(Hcount-1)));

% adjust threshold
if (MaximumT1WM-LowThreshold)>=(HighThreshold-MaximumT1WM)
    LowThreshold = MaximumT1WM - (HighThreshold-MaximumT1WM);
else 
    HighThreshold = MaximumT1WM + (MaximumT1WM-LowThreshold);
end

% save threshold
%Threshold = [LowThreshold, HighThreshold];

%DOG SPECIFIC Grady, Yong, Kevin
% Threshold = [672, 840];
% % YIJ 20170508: threshold test
%LowThreshold = 620;
%HighThreshold = 820;
Threshold = [LowThreshold, HighThreshold];

% make WM mask
ROIsMaskWM = ones(size(T1map));
ROIsMaskWM(find(T1map<LowThreshold))=0;
ROIsMaskWM(find(T1map>HighThreshold))=0;
ROIsMaskWM=ROIsMaskWM.*mask;

% % remove periphery of mask
% ROIsMaskWM = fCutEdgeofMask(ROIsMaskWM, cutNpixel);
ROIsMaskWM = (outerbwperim(ROIsMaskWM)+ ROIsMaskWM) - bwperim(outerbwperim(ROIsMaskWM)+ ROIsMaskWM);
ROIsMaskWM = ROIsMaskWM - bwperim(ROIsMaskWM);

% clear variables
clear T1map Tesla mask appoxHighT1WM tempthreshold dT1  temp Index nMaxT1WM b nIndex MaximumT1WM HighThreshold LowThreshold count
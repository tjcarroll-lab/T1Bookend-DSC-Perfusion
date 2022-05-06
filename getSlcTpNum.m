function [numslc,numtime] = getSlcTpNum(dscpath)
%getSlcTpNum Get number of slices and timepoints from dsc
%   Inputs:
%       - dscpath: path to dsc (..\path\to\ep2d_perf')
%
% Author: Yong Ik Jeong
% Date: 2019-03-13
% Changelog:
%   - 20200313 YIJ: Initial version

bpdir = dir(fullfile(dscpath,'*.dcm'));
timeNumList = [];
slcNumList = [];

for ii = 1:length(bpdir)
    tmphdr = dicominfo(fullfile(bpdir(ii).folder,bpdir(ii).name));
    %hdrlist{ii} = tmphdr;
    acqNumList(ii) = tmphdr.AcquisitionNumber;
    imgNumList(ii) = tmphdr.InstanceNumber;
    if isfield(tmphdr,'TemporalPositionIdentifier')
        timeNumList(ii) = tmphdr.TemporalPositionIdentifier;
    end
    slcNumList(ii) = tmphdr.SliceLocation;
end

slcuniq = unique(slcNumList);
timeuniq = unique(timeNumList);

if ~isempty(timeNumList)
    numslc = length(slcuniq);
    numtime = length(timeuniq);
    if numslc*numtime ~= length(bpdir)
        error('Number of timepoints and slices do not match total number of images.');
    end
else
    numslc = length(slcuniq);
    numtime = length(bpdir)/numslc;
    warning('Could not find number of time points from header, so estimating from total number of images.');
end

end


function [copypath] = sortPerf_ICAD_UC(bpcase)
%sortPerf_ICAD_UC Reorders perfusion series for compatibility with AutoqCBF
%   Auto qCBF code expects timeseries to be ordered by time first for each
%   slice (s0t0, s0t1, s0t2, ..., s1t0, s1t1, s1t2, ...)
%   Also use it for LLs to rename them in 1.dcm, 2.dcm, ... format
%
% Author: Yong Ik Jeong
% Date: 2020-02-27

%bp = 'D:\Users\CarrollLab\Desktop\CREST-H-UCLA sorted\Case 23';
bpdir = dir(fullfile(bpcase,'**/*.dcm'));
if isempty(bpdir)
    bpdir = dir(fullfile(bpcase,'IM*'));
end
%hdrlist = [];
tmpNumList = 1:length(bpdir);
acqNumList = [];
imgNumList = [];
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

if length(unique(acqNumList)) > 1
    sorttype = 'diffAcqNum';
else
    sorttype = '';
end

switch sorttype
    case 'diffAcqNum'
        
        sortInd = sortrows([slcNumList' imgNumList' acqNumList'  tmpNumList'],1,'descend');
        
        numtime = length(unique(acqNumList));
        numslc = length(bpdir)/numtime;
        if mod(length(bpdir),numtime) ~= 0
            error(sprintf('number of images does not match up with number of slices and timepoints (total images: %d, timepoints: %d)', length(bpdir), numtime));
        end

        sortInd = sortInd(:,end);
                     
    otherwise

        if isfield(tmphdr,'NumberOfTemporalPositions')
            numtime = tmphdr.NumberOfTemporalPositions;
        else
            numtime = length(unique(acqNumList));
        end
        
        sortInd = sortrows([slcNumList' timeNumList' imgNumList' acqNumList' tmpNumList'],1,'descend');

        numslc = length(bpdir)/numtime;
        if mod(length(bpdir),numtime) ~= 0
            error('number of images does not match up with number of slices and timepoints');
        end

        sortInd = sortInd(:,end);

end

bpdir = bpdir(sortInd);

for ii = 1:length(bpdir)
    
    copyfrom = fullfile(bpdir(ii).folder,bpdir(ii).name);

    copyto = fullfile(bpcase,[num2str(ii) '.dcm']);
    
    movefile(copyfrom,copyto);
end

copypath = fullfile(bpcase);
end


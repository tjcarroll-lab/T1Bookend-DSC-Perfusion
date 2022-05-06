function [copypath] = reversePerfSlcOrder(bpcase)
%reversePerfSlcOrder Some MYRIAD perfusion cases have slices ordered from
%top to bottom and conflicts with nifti conversion which sorts slices from
%bottom to top. So re-sort perfusion dicoms accordingly
%
% Author: Yong Ik Jeong
% Date: 2021-10-12

%bp = 'D:\Users\CarrollLab\Desktop\CREST-H-UCLA sorted\Case 23';
bpdir = dir(fullfile(bpcase,'**/*.dcm'));
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
        
        sortInd = sortrows([slcNumList' imgNumList' acqNumList'  tmpNumList'],[-1 2]);
        
        numtime = length(unique(acqNumList));
        numslc = length(bpdir)/numtime;
        if mod(length(bpdir),numtime) ~= 0
            error(sprintf('number of images does not match up with number of slices and timepoints (total images: %d, timepoints: %d)', length(bpdir), numtime));
        end
        %sortInd = [sortInd repmat([1:numslc]',numtime,1)];
        
        %sortInd = sortrows(sortInd(:,[end 1:end-1]));
        sortInd = sortInd(:,end);
                     
    otherwise

        if isfield(tmphdr,'NumberOfTemporalPositions')
            numtime = tmphdr.NumberOfTemporalPositions;
        else
            numtime = length(unique(acqNumList));
            if length(numtime) == 1
                error('Number of time is one?');
            end
        end
        
        sortInd = sortrows([slcNumList' timeNumList' imgNumList' acqNumList' tmpNumList'],[-1 2]);

        numslc = length(bpdir)/numtime;
        if mod(length(bpdir),numtime) ~= 0
            error('number of images does not match up with number of slices and timepoints');
        end
        %sortInd = [sortInd repmat([1:numslc]',numtime,1)];
        
        %sortInd = sortrows(sortInd(:,[end 1:end-1]));
        sortInd = sortInd(:,end);

end

bpdir = bpdir(sortInd);

perfdir = fullfile(bpcase,'..','ep2d_perf');
perfdirtmp = fullfile(bpcase,'..','ep2d_perf_temp');
movefile(perfdir,perfdirtmp);
mkdir(perfdir);
for ii = 1:length(bpdir)
    
    copyfrom = fullfile(perfdirtmp,bpdir(ii).name);
    copyto = fullfile(bpcase,'..','ep2d_perf',[num2str(ii) '.dcm']);
    
    copyfile(copyfrom,copyto);
end

copypath = fullfile(bpcase,'..','ep2d_perf');
end


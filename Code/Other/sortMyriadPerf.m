function [copypath] = sortMyriadPerf(bpcase)
%sortCrestHAcq Reorders perfusion series for compatibility with AutoqCBF
%   Auto qCBF code expects timeseries to be ordered by time first for each
%   slice (s0t0, s0t1, s0t2, ..., s1t0, s1t1, s1t2, ...)
%
% Author: Yong Ik Jeong
% Date: 2020-02-27

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
        
        sortInd = sortrows([slcNumList' imgNumList' acqNumList'  tmpNumList']);
        
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
        end
        
        sortInd = sortrows([slcNumList' timeNumList' imgNumList' acqNumList' tmpNumList']);

        numslc = length(bpdir)/numtime;
        if mod(length(bpdir),numtime) ~= 0
            error('number of images does not match up with number of slices and timepoints');
        end
        %sortInd = [sortInd repmat([1:numslc]',numtime,1)];
        
        %sortInd = sortrows(sortInd(:,[end 1:end-1]));
        sortInd = sortInd(:,end);

end

bpdir = bpdir(sortInd);

% 20210224 YIJ: commented out to use folder name instead
% try
%     if ~isfield(tmphdr,'SeriesDescription')
%         seriesname = tmphdr.ProtocolName;
%     else
%         seriesname = tmphdr.SeriesDescription;
%         if isempty(seriesname)
%             seriesname = tmphdr.ProtocolName;
%         end
%     end
% catch
%     seriesname = tmphdr.ProcedureCodeSequence.Item_1.CodeMeaning;
% end
% seriesname = formatSeries(seriesname);

for ii = 1:length(bpdir)
    
    copyfrom = fullfile(bpdir(ii).folder,bpdir(ii).name);
    
    % 20210224 YIJ: commented out to use folder name instead 
%     if ~exist(fullfile(bpcase,[seriesname ' FE_EPI']),'dir')
%         mkdir(fullfile(bpcase,[seriesname ' FE_EPI']));
%     end
%     copyto = fullfile(bpcase,[seriesname ' FE_EPI'],[seriesname ' FE_EPI_' num2str(ii) '.dcm']);
    if ~exist(fullfile(bpcase,'..','ep2d_perf'), 'dir')
        mkdir(fullfile(bpcase,'..','ep2d_perf'));
    end
    copyto = fullfile(bpcase,'..','ep2d_perf',[num2str(ii) '.dcm']);
    
    movefile(copyfrom,copyto);
end

copypath = fullfile(bpcase,'..','ep2d_perf');
end


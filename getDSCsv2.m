function [ dscs ] = getDSCsv2( bp, pickdate )
%getDSCsv2 load DSCs from dicoms (all timepoints per slice)
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2018-01-18

bp = checkFileSep(bp);
filelist = searchFolders(bp,'.dcm',0);

%pickid = {'03','04','05','06','07'};%,'10'};
%pickid = {'04','09','10','11'};
%pickid = {'10','11'};
%pickid = {'05'};
%pickid = {'03_test'};
pickid = {'04','05','06','12','21','23','26','30'};
if ~exist('pickdate','var')
    pickdate = {'d02'};
else
    if ~iscell(pickdate)
        pickdate = mat2cell(pickdate,1);
    end
end
pickresult = {'P001'};%,'P002'};
pickperf = {'ep2d_perf'};

dscs = struct([]);
count = 0;

for ii = 1:length(filelist)
    currpath = checkFileSep(filelist{ii});
    trailingpath = currpath(length(bp)+1:end);
    slashind = strfind(trailingpath,'\');
    if length(slashind) >= 4
        caseid = trailingpath(1:slashind(1)-1);
        casedate = trailingpath(slashind(1)+1:slashind(2)-1);
        caseresult = trailingpath(slashind(2)+1:slashind(3)-1);
        caseperf = trailingpath(slashind(3)+1:slashind(4)-1);
        
        if ~isempty(regexpi(caseid,'^COLLAT_')) && sum(strcmpi(caseid(8:end),pickid)) &&...
                sum(strcmpi(casedate(1:3),pickdate)) 
            if sum(strcmpi(caseresult,pickresult)) && sum(strcmpi(caseperf,pickperf))
                count = count + 1;
                dscs(count).id = caseid;
                dscs(count).date = casedate;
                
                tmphdr = dicominfo([currpath '\1.dcm']);
                ntp = tmphdr.NumberOfTemporalPositions;
                tmpdir = dir([currpath '\*dcm']);
                tmpdsc = [];
                tmpdsc_stack = [];
                for jj = 1:1:length(tmpdir)
                    tmpdsc = cat(3,tmpdsc,dicomread([currpath '\' num2str(jj) '.dcm']));
                    if mod(jj,ntp) == 0
                        tmpdsc_stack = cat(4,tmpdsc_stack,tmpdsc);
                        tmpdsc = [];
                    end
                end
                dscs(count).dsc_stack = tmpdsc_stack;
                dscs(count).header = tmphdr;
            end
        end
    end
end


end


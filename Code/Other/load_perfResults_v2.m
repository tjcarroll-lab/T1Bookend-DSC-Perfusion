function [ output ] = load_perfResults_v2( bp, pickdate )
%load_perfResults load perf results into struct variable
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2019-05-02

bp = checkFileSep(bp);
filelist = searchFolders(bp,'.mat',1);

%pickid = {'03','04','05','06','07','12'};
%pickid = {'04','05','06','12','21','23','26','30'};
pickid = {'07','13','14','15','17','19','22','24'};
%pickid = {'05','06','07'};
%pickid = {'04','06'};
%pickid = {'03_test'};
%pickid = {'04','05','06','07'};
if ~exist('pickdate','var')
    pickdate = {'d02'};
else
    if ~iscell(pickdate)
        pickdate = mat2cell(pickdate,1);
    end
end
    
pickresult = {'Result_MSwcf2'};
%pickperf = {'P001','P002'};
pickperf = {'P001GE_M.mat'};

output = struct([]);
count = 0;
count2 = 0;

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
                output(count).id = caseid;
                output(count).date = casedate;
                output(count).perf = caseperf;
                tmp = load(currpath(1:end-1));
                wmslice = tmp.ROIs.positions.n_slice_WM_DSC;
                tmpdsc = tmp.ROIs.dsc_stack;
                tmp.ROIs = rmfield(tmp.ROIs,'dsc_stack');
                tmp.ROIs.dsc_slice = tmpdsc{wmslice};
                tmp.ROIs = rmfield(tmp.ROIs,'sigmap');
                output(count).result = tmp;
                output(count).header = dicominfo([bp caseid '\' casedate '\' caseperf(1:4) '\ep2d_perf\1.dcm']);
            end
        end
    end
end


end


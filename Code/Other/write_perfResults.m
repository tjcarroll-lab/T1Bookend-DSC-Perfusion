function [ perfdata,perfstruct ] = write_perfResults( output, roi, write )
%write_perfResults write results to spreadsheet
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-07-17

if nargin < 3
    write = 0;
end
totalcases = length(output);
roiid = {roi.id};
count = 0;

for ii = 1:totalcases
    id = output(ii).id;
    if sum(strcmpi(id,roiid))
        image_names = output(ii).result.image_names;
        qcbf = output(ii).result.images{strmatch('qCBF_nSVD',image_names)};
%         if strcmpi(id,'COLLAT_07') && ~isempty(regexpi(output(ii).perf,'^P002'))
%             qcbf = qcbf(:,:,2:6);
%         end
        roi_stack = roi(strcmpi(id,roiid)).roi_stack;
        maxind = max(roi_stack(:));
        %perfdata = {};
        numbers = {};
        for jj = 1:maxind
            if sum(roi_stack == jj) == 0
                perfdata{jj,ii} = [];
            else
                perfdata{jj,ii} = mean(qcbf(roi_stack == jj));
            end
            numbers{jj,1} = jj;
        end
        
        if mod(ii,2)
            count = count + 1;
            cbfpre = perfdata(1:length(numbers),ii);
            perfstruct(count).casenum = str2double(regexpi(id,'\d\d','match'));
        else
            cbfpost = perfdata(1:length(numbers),ii);
            
            cond2 = cellfun(@isempty,cbfpre) | cellfun(@isempty,cbfpost);
            numbers(cond2) = [];
            cbfpre(cond2) = [];
            cbfpost(cond2) = [];
            
            cond1 = cellfun(@isnan,cbfpre) | cellfun(@isnan,cbfpost);
            cond3 = cell2mat(cellfun(@(x) x==0,cbfpre,'UniformOutput',false)) |...
                cell2mat(cellfun(@(x) x==0,cbfpost,'UniformOutput',false));
            numbers(cond1 | cond3) = [];
            cbfpre(cond1 | cond3) = [];
            cbfpost(cond1 | cond3) = [];
            
            perfstruct(count).cbfpre = [numbers, cbfpre];
            perfstruct(count).cbfpost = [numbers, cbfpost];
            cbfpre = {};
            cbfpost = {};
        end
    end
end

if write
    %xlswrite('..\DEBUG\perfdata_20170814_2.xlsx',perfdata);
    xlswrite(['..\DEBUG\' write],perfdata);
end

end

